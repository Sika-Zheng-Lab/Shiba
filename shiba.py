import argparse
import os
import sys
import logging
import yaml
import subprocess
# Configure logger
logger = logging.getLogger(__name__)
# Set version
VERSION = "v0.4.1"

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
    - Calculates PSI values and perform differential analysis
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
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")
    return parser.parse_args()

def check_file(file_path, description):
    if not os.path.isfile(file_path):
        logger.error(f"{description} does not exist: {file_path}")
        sys.exit(1)

def load_config(config_path):
    """
    Loads a YAML configuration file.

    Parameters:
    config_path (str): Path to the YAML configuration file.

    Returns:
    dict: The loaded configuration as a dictionary.
    """
    try:
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML configuration file: {e}")

def check_config(config, keys):
    for key in keys:
        if key not in config or not config[key]:
            logger.error(f"{key} is not set in the config file.")
            sys.exit(1)

def execute_command(command):
    logger.info(f"Executing: {' '.join(command)}")
    result = subprocess.run(command, shell=False)
    if result.returncode != 0:
        logger.error(f"Command failed: {' '.join(command)}")
        sys.exit(1)
    logger.info(f"Completed.")

def main():

    # Get arguments
    args = parse_args()
    # Set up logging
    logging.basicConfig(
        format = "[%(asctime)s] %(levelname)7s %(message)s",
        level = logging.DEBUG if args.verbose else logging.INFO
    )

    # Validate input and config
    logger.info("Running Shiba...")
    config_path = args.config
    processors = str(args.process)

    # Load config
    logger.info("Loading configuration...")
    config = load_config(config_path)

    # Check essential config keys
    check_config(config, ["workdir", "experiment_table", "gtf"])

    # Prepare output directory
    output_dir = config["workdir"]
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "Shiba.log")
    logging.FileHandler(log_file)

    # Get parent directory of this script
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Get command line of this script
    command_line = " ".join(sys.argv)

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
                "-p", processors,
                "-r", config['reference_group'],
                "-a", config['alternative_group'],
                "-f", str(config['fdr']),
                "-d", str(config['delta_psi']),
                "-m", str(config['minimum_reads']),
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
                "-r", config['reference_group'],
                "-a", config['alternative_group'],
                "-p", processors
            ]
        },
        {
            "name": "Step 6: pca.py",
            "command": [
                "python", os.path.join(script_dir, "src", "pca.py"),
                "--input-tpm", os.path.join(output_dir, "results", "expression", "TPM.txt"),
                "--input-psi", os.path.join(output_dir, "results", "splicing", "PSI_matrix_group.txt" if config['only_psi_group'] else "PSI_matrix_sample.txt"),
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

    # Get start step
    start_step = args.start_step
    if start_step > 0:
        logger.info(f"Starting from step {start_step}...")
        steps = steps[start_step-1:]

    # Execute steps
    for step in steps:
        logger.info(f"Executing {step['name']}...")
        command_to_run = step["command"] + ["-v"] if args.verbose else step["command"]
        logger.debug(command_to_run)
        execute_command(command_to_run)

    # Finish
    logger.info(f"Shiba finished! Results saved in {output_dir}")

if __name__ == "__main__":
    main()
