import argparse
import os
import sys
import logging
import subprocess
import yaml
from src.lib import general
import time
# Configure logger
logger = logging.getLogger(__name__)
# Set version
VERSION = "v0.5.2"

def parse_args():
    parser = argparse.ArgumentParser(
                description=f"""scShiba {VERSION} - Pipeline for identification of differential RNA splicing in single-cell RNA-seq data

Step 1: gtf2event.py
    - Converts GTF files to event format.
Step 2: sc2junc.py
    - Counts junction reads from STARsolo output files.
Step 3: scpsi.py
    - Calculates PSI values and perform differential analysis.""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("config", help="Config file in yaml format")
    parser.add_argument("-p", "--process", type=int, default=1, help="Number of processors to use (default: 1)")
    parser.add_argument("-s", "--start-step", type=int, default=0, help="Start the pipeline from the specified step (default: 0, run all steps)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
    return parser.parse_args()

def main():

	# Get arguments
	args = parse_args()

	# Set up logging
	logging.basicConfig(
		format = "[%(asctime)s] %(levelname)7s %(message)s",
		level = logging.DEBUG if args.verbose else logging.INFO
	)

	# Validate input and config
	logger.info(f"Running scShiba ({VERSION})")
	time.sleep(1)
	logger.debug(f"Arguments: {args}")
	# Get number of processors
	processors = str(args.process)

	# Load config
	logger.info("Loading configuration...")
	config_path = args.config
	config = general.load_config(config_path)

	# Check essential config keys
	missing_keys = general.check_config(config, ["workdir", "experiment_table", "gtf"])
	if missing_keys:
		logger.error(f"Missing required keys in configuration file: {', '.join(missing_keys)}")
		sys.exit(1)
	else:
		logger.info(f"workdir: {config['workdir']}")
		logger.info(f"experiment_table: {config['experiment_table']}")
		logger.info(f"gtf: {config['gtf']}")

	# Prepare output directory
	output_dir = config["workdir"]
	logger.debug("Making output directory...")
	os.makedirs(output_dir, exist_ok=True)

	# Get parent directory of this script
	logger.debug("Getting script directory...")
	script_dir = os.path.dirname(os.path.realpath(__file__))

	# Get command line of this script
	logger.debug("Getting command line...")
	command_line = " ".join(sys.argv)
	logger.debug(f"Command line: {command_line}")

	# Steps
	experiment_table = config["experiment_table"]
	gtf = config["gtf"]
	steps = [
		{
			"name": "Step 1: gtf2event.py",
			"command": [
				"python", os.path.join(script_dir, "src", "gtf2event.py"),
				"-i", gtf,
				"-o", os.path.join(output_dir, "events"),
				"-p", processors
			]
		},
		{
			"name": "Step 2: sc2junc.py",
			"command": [
				"python", os.path.join(script_dir, "src", "sc2junc.py"),
				"-i", experiment_table,
				"-o", os.path.join(output_dir, "junctions", "junctions.bed")
			]
		},
		{
			"name": "Step 3: scpsi.py",
			"command": [
				"python", os.path.join(script_dir, "src", "scpsi.py"),
				"-p", "1",
				"-r", config['reference_group'],
				"-a", config['alternative_group'],
				"-f", str(config['fdr']),
				"-d", str(config['delta_psi']),
				"-m", str(config['minimum_reads']),
				"--onlypsi" if config['only_psi'] else "",
				"--excel" if config['excel'] else "",
				os.path.join(output_dir, "junctions", "junctions.bed"),
				os.path.join(output_dir, "events"),
				os.path.join(output_dir, "results")
			]
		}
	]

	logger.info("Starting pipeline execution...")

	# Get start step
	start_step = args.start_step
	if start_step > 0:
		logger.info(f"Starting from step {start_step}...")
		steps = steps[start_step-1:]

	# Execute steps
	for step in steps:
		logger.info(f"Executing {step['name']}...")
		command_to_run = step["command"] + ["-v"] if args.verbose else step["command"]
		# Delete empty strings
		command_to_run = [x for x in command_to_run if x]
		logger.debug(command_to_run)
		returncode = general.execute_command(command_to_run)
		if returncode != 0:
			logger.error(f"Error executing {step['name']}. Exiting...")
			sys.exit(1)

	# Finish
	logger.info(f"scShiba finished! Results saved in {output_dir}")
	# Generate report
	report_name = "scShiba"
	general.generate_report(report_name, output_dir, VERSION, command_line, experiment_table)

if __name__ == "__main__":
    main()
