import argparse
import os
import sys
import shutil
import subprocess
import logging
import pandas as pd
import pysam
from lib import expression, general

# Configure logging
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Shared helper functions
# ---------------------------------------------------------------------------

def prepare_output_dir(output_path):
	output_dir = os.path.dirname(output_path)
	logs_dir = os.path.join(output_dir, "logs")
	tmp_dir = os.path.join(output_dir, "tmp")
	os.makedirs(logs_dir, exist_ok=True)
	os.makedirs(tmp_dir, exist_ok=True)
	return output_dir, logs_dir, tmp_dir

def create_saf_file(ri_event, tmp_dir):
	saf_file = os.path.join(tmp_dir, "RI.saf")
	logger.info(f"Generating SAF file: {saf_file}")
	saf_data = []
	with open(ri_event, "r") as ri_file:
		next(ri_file)  # Skip header
		for line in ri_file:
			intron, strand = line.strip().split("\t")[5:7]
			chrom = intron.split(":")[0]
			start = intron.split(":")[1].split("-")[0]
			start_plus_1 = str(int(start) + 1)
			end = intron.split(":")[1].split("-")[1]
			end_minus_1 = str(int(end) - 1)
			saf_data.append([f"{chrom}:{start}-{start_plus_1}", chrom, start, start_plus_1, strand])
			saf_data.append([f"{chrom}:{end_minus_1}-{end}", chrom, end_minus_1, end, strand])
	saf_df = pd.DataFrame(saf_data, columns=["GeneID", "Chr", "Start", "End", "Strand"])
	saf_df.drop_duplicates().to_csv(saf_file, sep="\t", index=False)
	return saf_file

def run_featurecounts_ri(bam, ri_saf, output, threads, long_read=False):
	"""Run featureCounts for RI (exon-intron) junction counting on a single BAM.

	Used by both the all-in-one pipeline and the 'ri' subcommand.
	"""
	# Check if BAM is paired-end
	paired_flag = expression.is_paired_end(bam)
	paired_option = ["-p", "-B"] if paired_flag else []
	# Check if long read mode is enabled
	longread_option = ["-L"] if long_read else []
	# Run featureCounts
	featurecounts_command = [
		"featureCounts", "-a", ri_saf, "-o", output, "-F",
		"SAF", "--fracOverlapFeature", "1.0", "-T", str(threads), "-O"
	] + paired_option + longread_option + [bam]
	# Delete empty strings
	featurecounts_command = list(filter(None, featurecounts_command))
	returncode = general.execute_command(featurecounts_command)
	if returncode != 0:
		logger.error("Error executing featureCounts. Exiting...")
		sys.exit(1)

def process_samples(experiment_file, strand, anchor, min_intron, max_intron, output_dir, logs_dir, tmp_dir, saf_file, processors):
	junc_files = []
	with open(experiment_file, "r") as experiment:
		for line in experiment:
			line = line.strip()
			if not line or line.startswith("sample"):
				continue
			# Parse the line
			try:
				sample, bam, _group, technology = line.split(maxsplit=3)
			except ValueError:
				sample, bam, _group = line.split(maxsplit=2)
				technology = "short"
			sample_dir = os.path.dirname(bam)
			bam_index = f"{bam}.bai"

			logger.info(f"Processing sample: {sample}")
			logger.debug(f"BAM file: {bam}")
			if technology.lower() == "long":
				logger.debug(f"{sample} will be processed as a long read sequencing experiment.")

			# Ensure BAM index exists
			if not os.path.isfile(bam_index):
				logger.error(f"BAM index file not found for sample : {bam}")
				logger.error("Please create an index file using 'samtools index' for all BAM files.")
				sys.exit(1)
			else:
				logger.debug(f"Found BAM index for {bam}")

			# Extract exon-exon junctions
			logger.info(f"Counting exon-exon junctions for sample {sample}...")
			exon_junc_file = os.path.join(tmp_dir, f"{sample}_exon-exon.junc")
			regtools_command = [
				"regtools",
				"junctions",
				"extract",
				"-s", strand,
				"-a", str(anchor),
				"-m", str(min_intron),
				"-M", str(max_intron),
				"-o", exon_junc_file,
				bam
			]
			logger.debug(f"Regtools command: {regtools_command}")
			return_code = general.execute_command(
				regtools_command, os.path.join(logs_dir, "regtools.log")
			)
			if return_code != 0:
				logger.error(f"Regtools failed for sample {sample}")
				sys.exit(1)
			junc_files.append((exon_junc_file, "exon-exon"))

			# Count exon-intron junctions using shared helper
			logger.info(f"Counting exon-intron junctions for sample {sample}...")
			exon_intron_file = os.path.join(tmp_dir, f"{sample}_exon-intron.junc")
			long_read = technology.lower() == "long"
			run_featurecounts_ri(bam, saf_file, exon_intron_file, processors, long_read=long_read)
			junc_files.append((exon_intron_file, "exon-intron"))

	return junc_files

# ---------------------------------------------------------------------------
# Junction merge logic (shared by all-in-one and 'merge' subcommand)
# ---------------------------------------------------------------------------

def merge_exonexon(filelist):
	"""Merge exon-exon junction files into a single DataFrame."""
	result_df = pd.DataFrame()
	for i in range(len(filelist)):
		logger.debug(filelist[i].split("/")[-1] + "...")
		junc_df = pd.read_csv(
			filelist[i],
			sep="\t",
			header=None,
			dtype=str
		)
		sample = filelist[i].split("/")[-1].rstrip("_exon-exon.junc")
		junc_df = junc_df.iloc[:, [0, 1, 2, 4, 10]]
		junc_df.columns = ["chr", "start", "end", "count", "block"]
		junc_df = junc_df.reset_index()
		junc_df.loc[(junc_df["chr"].str.isdecimal() == True) | (junc_df["chr"].str.len() <= 2), "chr"] = "chr" + junc_df["chr"]
		junc_df["count"] = junc_df["count"].astype("int32")
		junc_df["blockSize1"] = junc_df["block"].str.split(",", expand=True)[0].astype("int32")
		junc_df["start"] = junc_df["start"].astype("int32") + junc_df["blockSize1"]
		junc_df["blockSize2"] = junc_df["block"].str.split(",", expand=True)[1].astype("int32")
		junc_df["end"] = junc_df["end"].astype("int32") - junc_df["blockSize2"] + 1
		junc_df["ID"] = junc_df["chr"].astype(str) + ":" + junc_df["start"].astype(str) + "-" + junc_df["end"].astype(str)
		junc_df["sample"] = sample
		junc_df = junc_df[["ID", "sample", "count"]]
		# Group by ID and sample
		junc_df = junc_df.groupby(["ID", "sample"], as_index=False).sum()
		result_df = pd.concat([result_df, junc_df], axis=0) if not result_df.empty else junc_df

	result_df["count"] = result_df["count"].astype("int32")
	# Check if there are duplicated junctions
	if result_df.duplicated(subset=["ID", "sample"]).any():
		duplicates = result_df[result_df.duplicated(subset=["ID", "sample"], keep=False)]
		logger.debug(f"Duplicated junctions found: {duplicates}")
		logger.debug("Duplicated junctions occur when the same junction is detected in both strands.")
		logger.debug("The duplicated junctions will be summed and merged.")
		result_df = result_df.groupby(["ID", "sample"], as_index=False).sum()

	result_df = result_df.pivot(
		index="ID",
		columns="sample",
		values="count"
	).fillna(0).reset_index()
	result_df = result_df.rename(columns={"index": "ID"})

	result_df["chr"] = result_df["ID"].str.split(":", expand=True)[0]
	result_df["start"] = result_df["ID"].str.split(":", expand=True)[1].str.split("-", expand=True)[0].astype("int32")
	result_df["end"] = result_df["ID"].str.split(":", expand=True)[1].str.split("-", expand=True)[1].astype("int32")
	result_df["chr-start"] = result_df["chr"] + "-" + result_df["start"].astype(str)
	result_df["chr-end"] = result_df["chr"] + "-" + result_df["end"].astype(str)
	col = [i for i in result_df.columns if i not in ["chr", "start", "end", "ID", "chr-start", "chr-end", "mean"]]
	for j in col:
		result_df = result_df.astype({j: "int32"})
	col = ["chr", "start", "end", "ID"] + col
	result_df = result_df[col]

	return result_df

def merge_exonintron(filelist):
	"""Merge exon-intron junction files into a single DataFrame."""
	result_df = pd.DataFrame()
	for i in range(len(filelist)):
		logger.debug(filelist[i].split("/")[-1] + "...")
		junc_df = pd.read_csv(
			filelist[i],
			sep="\t",
			comment="#",
			dtype=str
		)
		sample = filelist[i].split("/")[-1].rstrip("_exon-intron.junc")
		junc_df = junc_df.iloc[:, [1, 2, 3, 0, 6]]
		junc_df.columns = ["chr", "start", "end", "ID", "count"]
		junc_df["sample"] = sample
		junc_df.loc[(junc_df["chr"].str.isdecimal() == True) | (junc_df["chr"].str.len() <= 2), "chr"] = "chr" + junc_df["chr"]
		result_df = pd.concat([result_df, junc_df], axis=0) if result_df is not None else junc_df

	result_df["count"] = result_df["count"].astype("int32")
	result_df = result_df.pivot(
		index=["chr", "start", "end", "ID"],
		columns="sample",
		values="count"
	).fillna(0).reset_index()
	col = [i for i in result_df.columns if i not in ["chr", "start", "end", "ID"]]
	for j in col:
		result_df = result_df.astype({j: "int32"})

	return result_df

def merge_and_save(exonexon_files, exonintron_files, output_file):
	"""Merge exon-exon and exon-intron junction files and save result."""
	# exon-exon junctions
	logger.info("Merge exon-exon junction count...")
	exon_exon_junc_df = merge_exonexon(exonexon_files)

	# exon-intron junctions
	logger.info("Merge exon-intron junction count...")
	exon_intron_junc_df = merge_exonintron(exonintron_files)

	# Combine and save results
	logger.info("Combine and save results...")
	result_df = pd.concat(
		[exon_exon_junc_df, exon_intron_junc_df]
	).sort_values(["chr", "start"])
	result_df.to_csv(
		output_file,
		sep="\t",
		index=False
	)

	junc_num = str(result_df.count()[0])
	logger.debug(f"Total number of junctions: {junc_num}")
	logger.info("Merge junctions completed")

def merge_junction_files(junc_files, output_file):
	"""Merge junction files from all-in-one pipeline (list of tuples)."""
	exon_exon_files = [j[0] for j in junc_files if j[1] == "exon-exon"]
	exon_intron_files = [j[0] for j in junc_files if j[1] == "exon-intron"]
	merge_and_save(exon_exon_files, exon_intron_files, output_file)

# ---------------------------------------------------------------------------
# Subcommand: ri — run featureCounts for a single BAM (RI junctions)
# (replaces bam2junc_RI_snakemake.py)
# ---------------------------------------------------------------------------

def cmd_ri(args):
	logger.info("Running featureCounts for RI junctions...")
	run_featurecounts_ri(args.bam, args.RI, args.junc, args.threads, long_read=args.long_read)
	logger.info("Done.")

# ---------------------------------------------------------------------------
# Subcommand: merge — merge junction files
# (replaces merge_junc_snakemake.py)
# ---------------------------------------------------------------------------

def cmd_merge(args):
	logger.info("Starting merge junctions")
	merge_and_save(args.exonexon, args.exonintron, args.output)

# ---------------------------------------------------------------------------
# Default: all-in-one pipeline (original bam2junc.py behavior)
# ---------------------------------------------------------------------------

def cmd_all(args):
	output_dir, logs_dir, tmp_dir = prepare_output_dir(args.output)
	saf_file = create_saf_file(args.ri_event, tmp_dir)
	logger.info("Extracting junctions from BAM files...")
	junc_files = process_samples(
		args.input, args.strand, args.anchor, args.min_intron, args.max_intron, output_dir, logs_dir, tmp_dir, saf_file, args.processors
	)
	logger.debug(junc_files)
	logger.info("Merging junction read counts...")
	merge_junction_files(junc_files, args.output)

	# Cleanup
	logger.debug("Cleaning up temporary files...")
	shutil.rmtree(tmp_dir)  # Temporary directory

	logger.info("Junction read counts processing completed!")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
	parser = argparse.ArgumentParser(
		description="Pipeline for processing junction read counts."
	)
	subparsers = parser.add_subparsers(dest="subcommand")

	# Default (all-in-one) arguments — used when no subcommand is given
	parser.add_argument("-i", "--input", help="Experiment table")
	parser.add_argument("-r", "--ri_event", help="Intron retention event file")
	parser.add_argument("-o", "--output", help="Output junction read counts file")
	parser.add_argument("-p", "--processors", type=int, default=1, help="Number of processors to use (default: 1)")
	parser.add_argument("-a", "--anchor", type=int, default=8, help="Minimum anchor length (default: 8)")
	parser.add_argument("-m", "--min_intron", type=int, default=70, help="Minimum intron size (default: 70)")
	parser.add_argument("-M", "--max_intron", type=int, default=500000, help="Maximum intron size (default: 500000)")
	parser.add_argument("-s", "--strand", default="XS", help="Strand specificity (default: XS)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

	# Subcommand: ri
	ri_parser = subparsers.add_parser("ri", help="Run featureCounts for RI (exon-intron) junctions on a single BAM")
	ri_parser.add_argument("-b", "--bam", type=str, required=True, help="Input BAM file")
	ri_parser.add_argument("-r", "--RI", type=str, required=True, help="Input RI SAF file")
	ri_parser.add_argument("-o", "--junc", type=str, required=True, help="Output junction file")
	ri_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
	ri_parser.add_argument("-l", "--long-read", action="store_true", help="Long read mode")
	ri_parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

	# Subcommand: merge
	merge_parser = subparsers.add_parser("merge", help="Merge exon-exon and exon-intron junction files")
	merge_parser.add_argument("--exonexon", type=str, nargs="+", required=True, help="Exon-exon junction files")
	merge_parser.add_argument("--exonintron", type=str, nargs="+", required=True, help="Exon-intron junction files")
	merge_parser.add_argument("--output", type=str, required=True, help="Output file")
	merge_parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

	return parser.parse_args()

def main():
	args = parse_args()
	# Set up logging
	logging.basicConfig(
		format="[%(asctime)s] %(levelname)7s %(message)s",
		level=logging.DEBUG if args.verbose else logging.INFO
	)
	logger.debug(args)

	if args.subcommand == "ri":
		cmd_ri(args)
	elif args.subcommand == "merge":
		cmd_merge(args)
	else:
		# Default: all-in-one pipeline (original behavior)
		logger.info("Processing junction read counts...")
		cmd_all(args)

if __name__ == "__main__":
	main()
