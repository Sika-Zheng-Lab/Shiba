import argparse
import os
import sys
import subprocess
import logging
import pandas as pd
from lib import expression, general

# Configure logging
logger = logging.getLogger(__name__)

def parse_args():
	parser = argparse.ArgumentParser(
		description="RNA expression analysis using featureCounts and DESeq2."
	)
	parser.add_argument("-i", "--input", required=True, help="Experiment table")
	parser.add_argument("-g", "--reference", required=True, help="Reference GTF file")
	parser.add_argument("-o", "--output", required=True, help="Output directory")
	parser.add_argument("-r", "--refgroup", default="NA", help="Reference group for differential expression analysis")
	parser.add_argument("-a", "--altgroup", default="NA", help="Alternative group for differential expression analysis")
	parser.add_argument("-p", "--processors", type=int, default=1, help="Number of processors to use (default: 1)")
	parser.add_argument("--excel", help = "Make result files in excel format", action = 'store_true')
	parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
	return parser.parse_args()

def prepare_output_dir(output_dir):
    os.makedirs(f"{output_dir}/logs", exist_ok=True)

def process_samples(experiment_file, reference_gtf, output_dir, processors):
	count_all_df = pd.DataFrame()
	with open(experiment_file, "r") as experiment:
		for line in experiment:
			line = line.strip()
			if not line or line.startswith("sample"):
				continue
			# Parse the line
			try:
				sample, bam_file, _group, technology = line.split(maxsplit=3)
			except ValueError:
				sample, bam_file, _group = line.split(maxsplit=2)
				technology = "short"
			bam_index = f"{bam_file}.bai"

			logger.info(f"Processing sample: {sample}")
			logger.debug(f"BAM file: {bam_file}")
			if technology.lower() == "long":
				logger.debug(f"{sample} will be processed as a long read sequencing experiment.")

			# Ensure BAM index exists
			if not os.path.isfile(bam_index):
				logger.error(f"BAM index file not found for sample : {sample}")
				logger.error("Please create an index file using 'samtools index' for all BAM files.")
				sys.exit(1)
			else:
				logger.debug(f"Found BAM index for {bam_file}")

			# Check if BAM is paired-end
			paired_flag = expression.is_paired_end(bam_file)
			paired_option = ["-p", "-B"] if paired_flag else []
			# Check if BAM is long-read
			longread_option = ["-L"] if technology.lower() == "long" else []

			# Run featureCounts
			counts_file = f"{output_dir}/{sample}_counts.txt"
			featurecounts_command = [
				"featureCounts", "-a", reference_gtf, "-o", counts_file, "-T", str(processors),
				"-t", "exon", "-g", "gene_id"
			] + paired_option + longread_option + [bam_file]
			# Delete empty strings
			featurecounts_command = list(filter(None, featurecounts_command))
			return_code = general.execute_command(featurecounts_command, f"{output_dir}/logs/featureCounts.log")
			if return_code != 0:
				logger.error(f"FeatureCounts failed for sample {sample}")
				sys.exit(1)

			# Simplify counts file
			logger.info(f"Simplifying counts file for sample {sample}")
			counts_file_df = pd.read_csv(counts_file, sep="\t", comment="#")
			counts_file_df = counts_file_df.iloc[:, [0, 5, 6]]
			counts_file_df.columns = ["GeneID", "Length", sample]
			count_all_df = pd.merge(count_all_df, counts_file_df[["GeneID", sample]], on = "GeneID") if not count_all_df.empty else counts_file_df

	count_all_df = count_all_df.rename(columns = {"GeneID": "gene_name"})
	count_all_df = count_all_df.sort_values("gene_name")
	return count_all_df

def run_deseq2(src_path, experiment_file, counts_file, refgroup, altgroup, output_dir):
	if refgroup != "NA" and altgroup != "NA":
		logger.info("Running differential expression analysis using DESeq2...")
		deseq2_command = [
			"Rscript", f"{src_path}/deseq2.R", experiment_file, counts_file, refgroup, altgroup, f"{output_dir}/DEG.txt"
		]
		return_code = general.execute_command(deseq2_command, f"{output_dir}/logs/DESeq2.log")
		if return_code != 0:
			logger.error("DESeq2 failed")
			sys.exit(1)

def main():

	# Parse arguments
	args = parse_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)
	logger.info("Starting RNA expression analysis...")
	logger.debug(args)

	# Import StyleFrame
	from styleframe import StyleFrame, Styler, utils

	# Prepare output directory
	prepare_output_dir(args.output)

	# Process samples and generate count files
	count_all_df = process_samples(args.input, args.reference, args.output, args.processors)
	# Save count files
	logger.info("Saving count files...")
	count_df = count_all_df.drop(columns = ["Length"])
	count_df.to_csv(f"{args.output}/counts.txt", sep = "\t", index = False)

	# Calculate TPM
	logger.info("Calculating TPM...")
	tpm_df = expression.ExpressionProcessor(count_all_df.copy()).TPM()
	# Save TPM
	tpm_df.to_csv(f"{args.output}/TPM.txt", sep = "\t", index = False)

	# Calculate CPM
	logger.info("Calculating CPM...")
	cpm_df = expression.ExpressionProcessor(count_all_df.copy()).CPM()
	# Save CPM
	cpm_df.to_csv(f"{args.output}/CPM.txt", sep = "\t", index = False)

	# Export to Excel
	if args.excel:
		logger.info("Exporting results to Excel...")
		# StyleFrame
		style = Styler(
			horizontal_alignment = utils.horizontal_alignments.left,
			border_type = utils.borders.default_grid,
			wrap_text = False
		)
		with StyleFrame.ExcelWriter(f"{args.output}/TPM_CPM.xlsx") as writer:
			tpm_sf = StyleFrame(tpm_df)
			tpm_sf.set_column_width(columns = tpm_df.columns, width = 20)
			tpm_sf.apply_column_style(cols_to_style = tpm_df.columns, styler_obj = style, style_header = True)
			tpm_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "TPM")
			cpm_sf = StyleFrame(cpm_df)
			cpm_sf.set_column_width(columns = cpm_df.columns, width = 20)
			cpm_sf.apply_column_style(cols_to_style = cpm_df.columns, styler_obj = style, style_header = True)
			cpm_sf.to_excel(writer, index = False, columns_and_rows_to_freeze = "B2", sheet_name = "CPM")

	# Run DESeq2 for differential expression analysis
	run_deseq2(os.path.dirname(__file__), args.input, f"{args.output}/counts.txt", args.refgroup, args.altgroup, args.output)

	# Cleanup
	for count_file in os.listdir(args.output):
		if count_file.endswith("_counts.txt") or count_file.endswith("_counts.txt.summary"):
			os.remove(os.path.join(args.output, count_file))

	logger.info("RNA expression analysis completed!")

if __name__ == "__main__":
	main()
