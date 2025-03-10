import argparse
import os
import sys
import subprocess
import logging
import pandas as pd
import pysam
from lib import expression, general

# Configure logging
logger = logging.getLogger(__name__)

def get_args():
	parser = argparse.ArgumentParser(
		description="Pipeline for processing junction read counts."
	)
	parser.add_argument("-i", "--input", required=True, help="Experiment table")
	parser.add_argument("-r", "--ri_event", required=True, help="Intron retention event file")
	parser.add_argument("-o", "--output", required=True, help="Output junction read counts file")
	parser.add_argument("-p", "--processors", type=int, default=1, help="Number of processors to use (default: 1)")
	parser.add_argument("-a", "--anchor", type=int, default=8, help="Minimum anchor length (default: 8)")
	parser.add_argument("-m", "--min_intron", type=int, default=70, help="Minimum intron size (default: 70)")
	parser.add_argument("-M", "--max_intron", type=int, default=500000, help="Maximum intron size (default: 500000)")
	parser.add_argument("-s", "--strand", default="XS", help="Strand specificity (default: XS)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
	return parser.parse_args()

def prepare_output_dir(output_path):
	output_dir = os.path.dirname(output_path)
	logs_dir = os.path.join(output_dir, "logs")
	os.makedirs(logs_dir, exist_ok=True)
	return output_dir, logs_dir

def create_saf_file(ri_event, output_dir):
	saf_file = os.path.join(output_dir, "RI.saf")
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

def process_samples(experiment_file, strand, anchor, min_intron, max_intron, output_dir, logs_dir, saf_file, processors):
	junc_files = []
	with open(experiment_file, "r") as experiment:
		for line in experiment:
			line = line.strip()
			if not line or line.startswith("sample"):
				continue
			sample, bam, _group = line.split(maxsplit=2)
			sample_dir = os.path.dirname(bam)
			bam_index = f"{bam}.bai"

			# Ensure BAM index exists
			if not os.path.isfile(bam_index):
				logger.error(f"BAM index file not found for sample : {bam}")
				logger.error("Please create an index file using 'samtools index' for all BAM files.")
				sys.exit(1)
			else:
				logger.debug(f"Found BAM index for {bam}")

			# Extract exon-exon junctions
			logger.info(f"Counting exon-exon junctions for sample {sample}...")
			exon_junc_file = os.path.join(sample_dir, f"{sample}_exon-exon.junc")
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

			# Count exon-intron junctions
			logger.info(f"Counting exon-intron junctions for sample {sample}...")
			exon_intron_file = os.path.join(sample_dir, f"{sample}_exon-intron.junc")
			# Check if BAM is paired-end
			paired_flag = expression.is_paired_end(bam)
			paired_option = ["-p"] if paired_flag else [""]
			featurecounts_command = [
				"featureCounts",
				"-a", saf_file,
				"-o", exon_intron_file,
				"-F", "SAF",
				"--fracOverlapFeature", "1.0",
				"-T", str(processors),
				"-O"
			] + paired_option + [bam]
			# Delete empty strings
			featurecounts_command = list(filter(None, featurecounts_command))
			return_code = general.execute_command(
				featurecounts_command, os.path.join(logs_dir, "featureCounts.log")
			)
			if return_code != 0:
				logger.error(f"FeatureCounts failed for sample {sample}")
				sys.exit(1)
			junc_files.append((exon_intron_file, "exon-intron"))

	return junc_files

def merge_junction_files(junc_files, output_file):
	def process_junction_files(files, junction_type):
		logger.info(f"Merging {junction_type} junction files...")
		result_df = []
		for file in files:
			sample_name = os.path.basename(file).rsplit("_", 1)[0]
			logger.debug(f"Processing sample: {sample_name}")

			df = pd.read_csv(
				file,
				sep="\t",
				comment="#" if junction_type == "exon-intron" else None,
				header=None,
				dtype={0: str}
			)

			if junction_type == "exon-exon":
				df = df.iloc[:, [0, 1, 2, 4, 10]]
				df.columns = ["chr", "start", "end", "count", "block"]
				df = df.reset_index()
				df["chr"] = df["chr"].apply(lambda x: f"chr{x}" if x.isdecimal() or len(x) <= 2 else x)
				df["blockSize1"] = df["block"].str.split(",", expand = True)[0].astype("int32")
				df["start"] = df["start"].astype("int32") + df["blockSize1"]
				df["blockSize2"] = df["block"].str.split(",", expand = True)[1].astype("int32")
				df["end"] = df["end"].astype("int32") - df["blockSize2"] + 1
				df["ID"] = df["chr"].astype(str) + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
			else:
				df = df.iloc[1:, [1, 2, 3, 0, 6]]
				df.columns = ["chr", "start", "end", "ID", "count"]

			df["chr"] = df["chr"].apply(lambda x: f"chr{x}" if x.isdecimal() or len(x) <= 2 else x)
			df["sample"] = sample_name
			result_df.append(df[["ID", "sample", "count"]] if junction_type == "exon-exon" else df)

		merged_df = pd.concat(result_df, ignore_index=True).drop_duplicates()
		merged_df["count"] = merged_df["count"].astype(int)
		return merged_df

	# Separate files by junction type
	exon_exon_files = [j[0] for j in junc_files if j[1] == "exon-exon"]
	exon_intron_files = [j[0] for j in junc_files if j[1] == "exon-intron"]

	# Process exon-exon junctions
	exon_exon_df = process_junction_files(exon_exon_files, "exon-exon")
	if exon_exon_df.duplicated(subset=["ID", "sample"]).any():
		duplicates = exon_exon_df[exon_exon_df.duplicated(subset=["ID", "sample"], keep=False)]
		logger.debug(f"Duplicated junctions found: {duplicates}")
		logger.debug("Duplicated junctions occur when the same junction is detected in both strands.")
		logger.debug("The duplicated junctions will be summed and merged.")
		# Group by ID and sample and sum counts
		exon_exon_df = exon_exon_df.groupby(["ID", "sample"], as_index=False).sum()

	exon_exon_df = exon_exon_df.pivot(index="ID", columns="sample", values="count").fillna(0).reset_index()
	exon_exon_df["chr"], exon_exon_df["start"], exon_exon_df["end"] = zip(*exon_exon_df["ID"].str.extract(r'([^:]+):(\d+)-(\d+)').values)
	exon_exon_df = exon_exon_df.astype({"start": int, "end": int})
	exon_exon_df = exon_exon_df[["chr", "start", "end", "ID"] + [col for col in exon_exon_df.columns if col not in ["chr", "start", "end", "ID"]]]

	# Process exon-intron junctions
	exon_intron_df = process_junction_files(exon_intron_files, "exon-intron")
	exon_intron_df = exon_intron_df.pivot(index=["chr", "start", "end", "ID"], columns="sample", values="count").fillna(0).reset_index()
	exon_intron_df = exon_intron_df.astype({col: int for col in exon_intron_df.columns if col not in ["chr", "start", "end", "ID"]})

	# Combine and save results
	final_df = pd.concat([exon_exon_df, exon_intron_df], ignore_index=True).sort_values(["chr", "start"])
	# Make sure values are all integers
	final_df = final_df.astype({col: int for col in final_df.columns if col not in ["chr", "start", "end", "ID"]})
	final_df.to_csv(output_file, sep="\t", index=False)
	logger.info(f"Junction read counts merged into {output_file}")

def main():

	# Parse arguments
	args = get_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)
	logger.info("Processing junction read counts...")
	logger.debug(args)

	output_dir, logs_dir = prepare_output_dir(args.output)
	saf_file = create_saf_file(args.ri_event, output_dir)
	logger.info("Extracting junctions from BAM files...")
	junc_files = process_samples(
		args.input, args.strand, args.anchor, args.min_intron, args.max_intron, output_dir, logs_dir, saf_file, args.processors
	)
	logger.debug(junc_files)
	logger.info("Merging junction read counts...")
	merge_junction_files(junc_files, args.output)

	# Cleanup
	logger.debug("Cleaning up temporary files...")
	os.remove(saf_file) # SAF file
	for junc_file, junc_type in junc_files:
		os.remove(junc_file) # Junction files
		if junc_type == "exon-intron":
			os.remove(junc_file + ".summary")
	logger.info("Junction read counts processing completed!")

if __name__ == "__main__":
	main()
