import argparse
import os
import sys
import logging
from lib import general

# Configure logging
logger = logging.getLogger(__name__)

def parse_args():
	parser = argparse.ArgumentParser(
		description="Pipeline for transcript assembly using StringTie2"
	)
	parser.add_argument("-i", "--input", required=True, help="Experiment table")
	parser.add_argument("-r", "--reference", required=True, help="Reference GTF file")
	parser.add_argument("-o", "--output", required=True, help="Assembled GTF file")
	parser.add_argument("-p", "--processors", type=int, default=1, help="Number of processors to use (default: 1)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
	return parser.parse_args()

def validate_output_file(output_file):
	if os.path.isdir(output_file):
		logger.error(f"Output path is a directory, not a file: {output_file}")
		raise IsADirectoryError(f"{output_file} is a directory, not a file.")
	if os.path.exists(output_file) and not os.access(output_file, os.W_OK):
		logger.error(f"Output file is not writable: {output_file}")
		raise PermissionError(f"{output_file} is not writable.")

def main():

	# Parse arguments
	args = parse_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)

	experiment_file = args.input
	reference_gtf = args.reference
	assembled_gtf = args.output
	num_processors = args.processors

	logger.info("Starting transcript assembly...")
	logger.debug(args)

	# Validate output file
	validate_output_file(assembled_gtf)

	# Prepare output directory
	output_dir = os.path.dirname(assembled_gtf)
	os.makedirs(output_dir, exist_ok=True)

	gtf_list = []
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
			sample_dir = os.path.dirname(bam_file)
			bam_index = f"{bam_file}.bai"

			logger.info(f"Processing sample: {sample}")
			logger.debug(f"BAM file: {bam_file}")
			if technology.lower() == "long":
				logger.debug(f"{sample} will be processed as a long read sequencing experiment.")

			# Check if BAM index exists
			if not os.path.isfile(bam_index):
				logger.error(f"BAM index file not found for {bam_file}")
				logger.error(f"Please create an index file using 'samtools index' for all BAM files.")
				sys.exit(1)
			else:
				logger.debug(f"Found BAM index for {bam_file}")

			# Run StringTie2 for assembly
			sample_gtf = os.path.join(os.path.dirname(assembled_gtf), f"{sample}.assembled.gtf")
			longread_flag = ["-L"] if technology.lower() == "long" else []
			stringtie_command = [
				"stringtie",
				"-p", str(num_processors),
				"-G", reference_gtf,
				"-o", sample_gtf
			] + longread_flag + [bam_file]
			# Add -L option if long read experiment
			logger.debug(f"StringTie2 command: {stringtie_command}")
			return_code = general.execute_command(stringtie_command)
			if return_code != 0:
				logger.error(f"StringTie2 failed for sample {sample}")
				sys.exit(1)

			# Add to GTF list
			gtf_list.append(sample_gtf)

	# Merge GTF files
	logger.info("Merging GTF files...")
	merge_command = [
		"stringtie",
		"--merge",
		"-p", str(num_processors),
		"-G", reference_gtf,
		"-o", assembled_gtf
	] + gtf_list
	logger.debug(f"StringTie2 merge command: {merge_command}")
	return_code = general.execute_command(merge_command)
	if return_code != 0:
		logger.error("StringTie2 merge failed")
		sys.exit(1)

	# Cleanup assembled GTF files
	logger.info("Cleaning up intermediate GTF files...")
	for gtf_file in gtf_list:
		os.remove(gtf_file)

	# Done
	logger.info("Transcript assembly completed successfully!")

if __name__ == "__main__":

	main()
