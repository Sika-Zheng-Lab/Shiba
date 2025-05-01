import argparse
import os
import sys
import pysam
from lib import expression, general
import logging
# Configure logging
logger = logging.getLogger(__name__)

def parse_args():
	## Get arguments from command line

	parser = argparse.ArgumentParser(
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		description = "bam2junc_RI_snakemake.py"
	)
	parser.add_argument('-b', '--bam', type=str, help='Input bam file')
	parser.add_argument('-r', '--RI', type=str, help='Input RI file')
	parser.add_argument('-o', '--junc', type=str, help='Output junction file')
	parser.add_argument('-t', '--threads', type=int, help='Number of threads')
	parser.add_argument('-l', '--long-read', action='store_true', help='Long read mode')
	parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
	args = parser.parse_args()
	return args

def bam2junc(bam, RI, output, threads, long_read=False):

	# Check if BAM is paired-end
	paired_flag = expression.is_paired_end(bam)
	paired_option = ["-p", "-B"] if paired_flag else [""]

	# Check if long read mode is enabled
	longread_option = ["-L"] if long_read else [""]

	# Run featureCounts
	featurecounts_command = [
		"featureCounts", "-a", RI, "-o", output, "-F",
		"SAF", "--fracOverlapFeature", "1.0", "-T", str(threads), "-O"
	] + paired_option + longread_option + [bam]
	# Delete empty strings
	featurecounts_command = list(filter(None, featurecounts_command))
	returncode = general.execute_command(featurecounts_command)
	if returncode != 0:
		logger.error(f"Error executing the command. Exiting...")
		sys.exit(1)

def main():

	# Parse arguments
	args = parse_args()
	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)
	logger.debug(args)

	# Run featureCounts
	logger.info("Running featureCounts...")
	bam2junc(args.bam, args.RI, args.junc, args.threads, args.long_read)

	# Finish
	logger.info("Done.")

if __name__ == "__main__":

	main()
