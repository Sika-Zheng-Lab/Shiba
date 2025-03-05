# Modules used in shiba.py and scshiba.py
import os
import subprocess
import sys
import yaml
import logging
import datetime

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
    missing_keys = [key for key in keys if key not in config or not config[key]]
    return missing_keys

def execute_command(command, log_file=None):
    if log_file:
        with open(log_file, "a") as log:
            result = subprocess.run(command, shell=False, stdout=log, stderr=log)
    else:
        result = subprocess.run(command, shell=False)
    return result.returncode

def generate_report(name, output_dir, version, command_line, experiment_table):
    report_path = os.path.join(output_dir, "report.txt")
    with open(report_path, "w") as report_file:
        report_file.write(f"Pipeline: {name}\n")
        report_file.write(f"Version: {version}\n")
        report_file.write(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        report_file.write(f"Command Line: {command_line}\n")
        report_file.write(f"Experiment Table: {experiment_table}\n")
    logging.info(f"Report generated at {report_path}")
    return None
