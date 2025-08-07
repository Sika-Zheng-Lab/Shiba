#!/usr/bin/env python3

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
VERSION = "v0.7.0"

def parse_args():
    parser = argparse.ArgumentParser(
                description=f"""Shiba {VERSION} - Pipeline for identification of differential RNA splicing

Step 1: bam2gtf.py
    - Assembles transcript structures based on mapped reads using StringTie2.
Step 2: gtf2event.py
    - Converts GTF files to event format.
Step 3: bam2junc.py
    - Extracts junction reads from BAM files.
Step 4: psi.py
    - Calculates PSI values and performs differential analysis.
Step 5: expression.py
    - Analyzes gene expression.
Step 6: pca.py
    - Performs PCA.
Step 7: plots.py
    - Generates plots from results.""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("config", help="Config file in yaml format")
    parser.add_argument("-p", "--process", type=int, default=1, help="Number of processors to use (default: 1)")
    parser.add_argument("-s", "--start-step", type=int, default=0, help="Start the pipeline from the specified step (default: 0, run all steps)")
    parser.add_argument("--mame", action="store_true", help="Execute MameShiba, a lightweight version of Shiba, for only splicing analysis. Steps 5-7 will be skipped.")
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

    # Set up report name
    report_name = "MameShiba" if args.mame else "Shiba"

    # Validate input and config
    logger.info(f"Running {report_name} ({VERSION})")
    if args.mame:
        logger.info("Only splicing analysis will be performed.")
    time.sleep(1)
    logger.debug(f"Arguments: {args}")
    # Get number of processors
    processors = str(args.process)

    # Load config
    logger.info("Loading configuration...")
    config_path = args.config
    config = general.load_config(config_path)
    if config['excel'] and args.mame:
        logger.warning("Excel output is not available in MameShiba mode. Disabling excel option.")
        config['excel'] = False

    # Check essential config keys
    missing_keys = general.check_config(config, ["workdir", "experiment_table", "gtf"])
    if missing_keys:
        logger.error(f"Missing required keys in configuration file: {', '.join(missing_keys)}")
        sys.exit(1)
    else:
        logger.info(f"workdir: {config['workdir']}")
        logger.info(f"experiment_table: {config['experiment_table']}")
        logger.info(f"gtf: {config['gtf']}")

    # Check optional config keys
    only_psi = config.get('only_psi', False)
    only_psi_group = config.get('only_psi_group', False)

    # Validate sample and group sizes in the experiment table
    sample_count = general.check_samplesize(config["experiment_table"])
    group_count = general.check_groupsize(config["experiment_table"])

    if sample_count <= 0:
        logger.error("No samples found in experiment table. Exiting...")
        sys.exit(1)
    elif sample_count == 1:
        logger.info("Only one sample found in experiment table.")
        only_psi = True
        logger.info("Differential analysis mode disabled. Set only_psi to True.")

    if group_count <= 0:
        logger.error("No groups found in experiment table. Exiting...")
        sys.exit(1)
    elif group_count == 1:
        logger.info("Only one group found in experiment table.")
        if not (only_psi or only_psi_group):
            logger.info("Differential analysis mode disabled.")
            only_psi = True

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
            "name": "Step 1: bam2gtf.py",
            "command": [
                "python", os.path.join(script_dir, "src", "bam2gtf.py"),
                "-i", experiment_table,
                "-r", gtf,
                "-o", os.path.join(output_dir, "annotation", "assembled_annotation.gtf"),
                "-p", processors
            ]
        },
        {
            "name": "Step 2: gtf2event.py",
            "command": [
                "python", os.path.join(script_dir, "src", "gtf2event.py"),
                "-i", os.path.join(output_dir, "annotation", "assembled_annotation.gtf") if config['unannotated'] else gtf,
                "-r" if config['unannotated'] else "",
                gtf if config['unannotated'] else "",
                "-o", os.path.join(output_dir, "events"),
                "-p", processors
            ]
        },
        {
            "name": "Step 3: bam2junc.py",
            "command": [
                "python", os.path.join(script_dir, "src", "bam2junc.py"),
                "-i", experiment_table,
                "-r", os.path.join(output_dir, "events", "EVENT_RI.txt"),
                "-o", os.path.join(output_dir, "junctions", "junctions.bed"),
                "-p", processors,
                "-a", str(config['minimum_anchor_length']),
                "-m", str(config['minimum_intron_length']),
                "-M", str(config['maximum_intron_length']),
                "-s", config['strand']
            ]
        },
        {
            "name": "Step 4: psi.py",
            "command": [
                "python", os.path.join(script_dir, "src", "psi.py"),
                "-g", experiment_table,
                "-p", "1",
                "-r", config['reference_group'],
                "-a", config['alternative_group'],
                "-f", str(config['fdr']),
                "-d", str(config['delta_psi']),
                "-m", str(config['minimum_reads']),
                "-i" if config['individual_psi'] else "",
                "-t" if config['ttest'] else "",
                "--excel" if config['excel'] else "",
                "--onlypsi" if only_psi else "",
                "--onlypsi-group" if only_psi_group else "",
                os.path.join(output_dir, "junctions", "junctions.bed"),
                os.path.join(output_dir, "events"),
                os.path.join(output_dir, "results", "splicing")
            ]
        },
        {
            "name": "Step 5: expression.py",
            "command": [
                "python", os.path.join(script_dir, "src", "expression.py"),
                "-i", experiment_table,
                "-g", gtf,
                "-o", os.path.join(output_dir, "results", "expression"),
                "" if only_psi or only_psi_group else "-r",
                "" if only_psi or only_psi_group else config['reference_group'],
                "" if only_psi or only_psi_group else "-a",
                "" if only_psi or only_psi_group else config['alternative_group'],
                "--excel" if config['excel'] else "",
                "-p", processors
            ]
        },
        {
            "name": "Step 6: pca.py",
            "command": [
                "python", os.path.join(script_dir, "src", "pca.py"),
                "--input-tpm", os.path.join(output_dir, "results", "expression", "TPM.txt"),
                "--input-psi", os.path.join(output_dir, "results", "splicing", "PSI_matrix_group.txt" if only_psi_group else "PSI_matrix_sample.txt"),
                "-g", "3000",
                "-o", os.path.join(output_dir, "results", "pca"),
            ]
        },
        {
            "name": "Step 7: plots.py",
            "command": [
                "python", os.path.join(script_dir, "src", "plots.py"),
                "-i", os.path.join(output_dir, "results"),
                "-e", experiment_table,
                "-s", command_line,
                "-o", os.path.join(output_dir, "plots")
            ]
        }
    ]

    logger.info("Starting pipeline execution...")

    # Skip after Step 4 if MameShiba is enabled and disable excel option in config
    if args.mame:
        steps = steps[:4]

    # Get start step
    start_step = args.start_step
    if start_step > 0:
        logger.info(f"Starting from step {start_step}...")
        steps = steps[start_step-1:]

    # Skip plots.py step when only_psi or only_psi_group is True
    if not args.mame and (only_psi or only_psi_group):
        steps = steps[:-1]

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
            logger.error("You could try to run the command with --verbose to catch the more detailed error message.")
            logger.error("If you think this is a bug, please report it at https://github.com/Sika-Zheng-Lab/Shiba/issues")
            sys.exit(1)

    # Finish
    logger.info(f"{report_name} finished! Results saved in {output_dir}")
    # Generate report
    general.generate_report(report_name, output_dir, VERSION, command_line, experiment_table)

if __name__ == "__main__":
    main()
