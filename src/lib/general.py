# Modules used in shiba.py and scshiba.py
import os
import subprocess
import sys
import yaml
import json
import logging
import datetime
logger = logging.getLogger(__name__)

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

def generate_report(name, output_dir, version, command_line, experiment_table, start_time=None):
    """
    Generate a JSON report with pipeline execution information.
    
    Parameters:
    name (str): Name of the pipeline
    output_dir (str): Directory to save the report
    version (str): Version of the tool
    command_line (str): Command line used to run the pipeline
    experiment_table (str): Path to the experiment table file
    start_time (datetime): Start time of the pipeline execution (optional)
    
    Returns:
    str: Path to the generated report.json file
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get end time
    end_time = datetime.datetime.now()
    
    # Calculate duration
    if start_time is None:
        logger.warning("Start time not provided. Duration will be set to 0.")
        duration_seconds = 0
    else:
        duration_seconds = (end_time - start_time).total_seconds()
        
    # Create report data structure
    report_data = {
        "tool": {
            "name": name,
            "version": version
        },
        "run": {
            "command": command_line,
            "start_time": start_time.isoformat(timespec="microseconds") if start_time else 'NA',
            "end_time": end_time.isoformat(timespec="microseconds") if end_time else 'NA',
            "duration_seconds": round(duration_seconds, 3) if duration_seconds != 0 else 'NA'
        },
        "experiment": {
            "table_path": experiment_table
        }
    }
    
    # Write JSON report
    report_path = os.path.join(output_dir, "report.json")
    try:
        with open(report_path, "w", encoding="utf-8") as report_file:
            json.dump(report_data, report_file, ensure_ascii=False, indent=2)
    except IOError as e:
        logger.error(f"Failed to write report file: {e}")
        raise
    
    logger.info(f"Report (JSON) generated at {report_path}")
    return report_path

def check_samplesize(experiment_table):
    # Check if experiment table exists
    try:
        with open(experiment_table, "r") as table:
            sample_count = sum(1 for _ in table) - 1
    except FileNotFoundError:
        raise FileNotFoundError(f"Experiment table not found: {experiment_table}")
    except Exception as e:
        raise ValueError(f"Error reading experiment table: {e}")
    # Check if there are samples in the experiment table
    if sample_count <= 0:
        raise ValueError('No samples found in experiment table')
    else:
        logger.debug(f'{sample_count} samples found in experiment table')
    return sample_count

def check_groupsize(experiment_table):
    # Check if experiment table exists
    try:
        with open(experiment_table, "r") as table:
            groups = {}
            for i, line in enumerate(table):
                if i == 0:  # Skip header
                    continue
                columns = line.strip().split("\t")
                if len(columns) < 3:
                    raise ValueError(f"Invalid format in experiment table at line {i + 1}")
                group = columns[2]
                groups[group] = groups.get(group, 0) + 1
    except FileNotFoundError:
        raise FileNotFoundError(f"Experiment table not found: {experiment_table}")
    except Exception as e:
        raise ValueError(f"Error reading experiment table: {e}")
    # Check if there are groups in the experiment table
    if not groups:
        raise ValueError('No groups found in experiment table')
    else:
        logger.debug(f'{len(groups.keys())} groups found in experiment table')
        for group, size in groups.items():
            logger.debug(f'Group "{group}" has {size} samples')
    return len(groups.keys())
