# Modules used in shiba.py and scshiba.py
import argparse
import os
import subprocess
import sys
import yaml
import json
import logging
import datetime
logger = logging.getLogger(__name__)


def str2bool(v):
    """Convert string to boolean for argparse compatibility.

    Accepts both bare flags (via nargs='?' + const=True) and explicit
    string values passed by Snakemake config interpolation.

    Accepted truthy values: "True", "true", "1", "yes"
    Accepted falsy values:  "False", "false", "0", "no"

    Usage in argparse:
        parser.add_argument("--flag", type=str2bool, nargs="?", const=True, default=False)
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("true", "1", "yes"):
        return True
    elif v.lower() in ("false", "0", "no"):
        return False
    else:
        raise argparse.ArgumentTypeError(f"Boolean value expected, got '{v}'")

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
    start_time (datetime or str): Start time of the pipeline execution (optional).
        Accepts a datetime object or an ISO 8601 format string.
    
    Returns:
    str: Path to the generated report.json file
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse start_time if given as ISO format string
    if isinstance(start_time, str):
        try:
            start_time = datetime.datetime.fromisoformat(start_time)
        except (ValueError, TypeError):
            logger.warning(f"Invalid start_time string: '{start_time}'. Treating as None.")
            start_time = None
    
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

# ==============================================================================
# Configuration validation functions
# ==============================================================================

def validate_file_exists(path, label):
    """
    Validates that a file exists at the given path.

    Parameters:
    path (str): The file path to check.
    label (str): A human-readable label for error messages (e.g. "GTF file").

    Returns:
    list: A list of error message strings (empty if the file exists).
    """
    errors = []
    if not os.path.isfile(path):
        errors.append(f'{label} not found: {path}')
    return errors

def validate_experiment_table_columns(experiment_table, required_columns):
    """
    Validates that the experiment table contains the required columns.

    Parameters:
    experiment_table (str): Path to the experiment table TSV file.
    required_columns (list): List of required column names.

    Returns:
    list: A list of error message strings (empty if all columns are present).
    """
    errors = []
    try:
        with open(experiment_table, "r") as f:
            header_line = f.readline().strip()
            if not header_line:
                errors.append(f'Experiment table is empty: {experiment_table}')
                return errors
            columns = header_line.split("\t")
            missing = [col for col in required_columns if col not in columns]
            if missing:
                errors.append(
                    f'Experiment table {experiment_table} is missing required column(s): {", ".join(missing)}. '
                    f'Found columns: {", ".join(columns)}'
                )
    except FileNotFoundError:
        errors.append(f'Experiment table not found: {experiment_table}')
    return errors

def validate_groups_bulk(experiment_table, reference_group, alternative_group):
    """
    Validates that reference_group and alternative_group exist in the 'group'
    column of a bulk experiment table, and that they are not identical.

    Parameters:
    experiment_table (str): Path to the bulk experiment table TSV file.
    reference_group (str): The reference group name from config.
    alternative_group (str): The alternative group name from config.

    Returns:
    list: A list of error message strings (empty if valid).
    """
    errors = []
    # Check reference != alternative
    if reference_group == alternative_group:
        errors.append(
            f'reference_group and alternative_group must be different, '
            f'but both are "{reference_group}"'
        )
    # Read groups from experiment table
    try:
        available_groups = set()
        with open(experiment_table, "r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                columns = line.strip().split("\t")
                if len(columns) >= 3:
                    available_groups.add(columns[2])
        if reference_group not in available_groups:
            errors.append(
                f'reference_group "{reference_group}" not found in experiment table group column. '
                f'Available groups: {", ".join(sorted(available_groups))}'
            )
        if alternative_group not in available_groups:
            errors.append(
                f'alternative_group "{alternative_group}" not found in experiment table group column. '
                f'Available groups: {", ".join(sorted(available_groups))}'
            )
    except FileNotFoundError:
        errors.append(f'Experiment table not found: {experiment_table}')
    return errors

def validate_groups_sc(experiment_table, reference_group, alternative_group):
    """
    Validates that reference_group and alternative_group exist in the 'group'
    column of the barcode TSV files referenced by the scShiba experiment table,
    and that they are not identical.

    Parameters:
    experiment_table (str): Path to the scShiba experiment table TSV file.
    reference_group (str): The reference group name from config.
    alternative_group (str): The alternative group name from config.

    Returns:
    list: A list of error message strings (empty if valid).
    """
    errors = []
    # Check reference != alternative
    if reference_group == alternative_group:
        errors.append(
            f'reference_group and alternative_group must be different, '
            f'but both are "{reference_group}"'
        )
    # Read barcode file paths from experiment table
    try:
        barcode_files = []
        with open(experiment_table, "r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                columns = line.strip().split("\t")
                if len(columns) >= 1 and columns[0]:
                    barcode_files.append(columns[0])
        # Read groups from each barcode file
        available_groups = set()
        for barcode_file in barcode_files:
            if not os.path.isfile(barcode_file):
                errors.append(f'Barcode file not found: {barcode_file}')
                continue
            with open(barcode_file, "r") as f:
                for j, line in enumerate(f):
                    if j == 0:
                        continue
                    columns = line.strip().split("\t")
                    if len(columns) >= 2:
                        available_groups.add(columns[1])
        if available_groups:
            if reference_group not in available_groups:
                errors.append(
                    f'reference_group "{reference_group}" not found in barcode file group column. '
                    f'Available groups: {", ".join(sorted(available_groups))}'
                )
            if alternative_group not in available_groups:
                errors.append(
                    f'alternative_group "{alternative_group}" not found in barcode file group column. '
                    f'Available groups: {", ".join(sorted(available_groups))}'
                )
    except FileNotFoundError:
        errors.append(f'Experiment table not found: {experiment_table}')
    return errors

def normalize_regtools_strand(strand_value):
    """Normalize Shiba strand configuration to regtools-compatible values.

    regtools >= 1.0.0 expects XS, RF, or FR instead of the legacy numeric values
    0, 1, and 2. This helper keeps the older config values accepted while
    translating them to the newer spellings before invoking regtools.
    """
    if strand_value is None:
        return "XS"

    strand = str(strand_value).strip()
    legacy_mapping = {
        "0": "XS",
        "1": "RF",
        "2": "FR",
    }
    if strand in legacy_mapping:
        return legacy_mapping[strand]
    if strand in {"XS", "RF", "FR"}:
        return strand
    raise ValueError(f"Unsupported strand value: {strand_value!r}")


def validate_config_types(config, mode="bulk"):
    """
    Validates config parameter types and value ranges.

    Parameters:
    config (dict): The loaded configuration dictionary.
    mode (str): "bulk" for Shiba/SnakeShiba, "sc" for scShiba/SnakeScShiba.

    Returns:
    list: A list of error message strings (empty if all valid).
    """
    errors = []

    # fdr: float, 0 < fdr <= 1
    if 'fdr' in config:
        try:
            fdr = float(config['fdr'])
            if not (0 < fdr <= 1):
                errors.append(f'fdr must be between 0 (exclusive) and 1 (inclusive), got {config["fdr"]}')
        except (TypeError, ValueError):
            errors.append(f'fdr must be a number, got "{config["fdr"]}"')

    # delta_psi: float, 0 <= delta_psi <= 1
    if 'delta_psi' in config:
        try:
            delta_psi = float(config['delta_psi'])
            if not (0 <= delta_psi <= 1):
                errors.append(f'delta_psi must be between 0 and 1, got {config["delta_psi"]}')
        except (TypeError, ValueError):
            errors.append(f'delta_psi must be a number, got "{config["delta_psi"]}"')

    # minimum_reads: int, > 0
    if 'minimum_reads' in config:
        try:
            minimum_reads = int(config['minimum_reads'])
            if minimum_reads <= 0:
                errors.append(f'minimum_reads must be a positive integer, got {config["minimum_reads"]}')
        except (TypeError, ValueError):
            errors.append(f'minimum_reads must be an integer, got "{config["minimum_reads"]}"')

    # Bulk-specific validations
    if mode == "bulk":
        # minimum_anchor_length: int, > 0
        if 'minimum_anchor_length' in config:
            try:
                val = int(config['minimum_anchor_length'])
                if val <= 0:
                    errors.append(f'minimum_anchor_length must be a positive integer, got {config["minimum_anchor_length"]}')
            except (TypeError, ValueError):
                errors.append(f'minimum_anchor_length must be an integer, got "{config["minimum_anchor_length"]}"')

        # minimum_intron_length: int, > 0
        if 'minimum_intron_length' in config:
            try:
                val = int(config['minimum_intron_length'])
                if val <= 0:
                    errors.append(f'minimum_intron_length must be a positive integer, got {config["minimum_intron_length"]}')
            except (TypeError, ValueError):
                errors.append(f'minimum_intron_length must be an integer, got "{config["minimum_intron_length"]}"')

        # maximum_intron_length: int, > 0 and > minimum_intron_length
        if 'maximum_intron_length' in config:
            try:
                max_val = int(config['maximum_intron_length'])
                if max_val <= 0:
                    errors.append(f'maximum_intron_length must be a positive integer, got {config["maximum_intron_length"]}')
                elif 'minimum_intron_length' in config:
                    try:
                        min_val = int(config['minimum_intron_length'])
                        if max_val <= min_val:
                            errors.append(
                                f'maximum_intron_length ({max_val}) must be greater than '
                                f'minimum_intron_length ({min_val})'
                            )
                    except (TypeError, ValueError):
                        pass  # Already reported above
            except (TypeError, ValueError):
                errors.append(f'maximum_intron_length must be an integer, got "{config["maximum_intron_length"]}"')

        # strand: one of XS, RF, FR, 0, 1, 2 (legacy values accepted for compatibility)
        if 'strand' in config:
            valid_strands = {"XS", "RF", "FR", "0", "1", "2"}
            if str(config['strand']) not in valid_strands:
                errors.append(
                    f'strand must be one of {", ".join(sorted(valid_strands))}, '
                    f'got "{config["strand"]}"'
                )

    return errors

def apply_config_defaults(config):
    """
    Applies default values for optional configuration keys that may be missing.
    Mutates the config dict in-place and logs warnings for any keys that were set
    to their default values.

    Parameters:
    config (dict): The loaded configuration dictionary.
    """
    if 'beta_binomial' not in config:
        config['beta_binomial'] = False
        logger.warning("'beta_binomial' is not specified in the configuration file. Defaulting to False.")

def validate_config(config, mode="bulk"):
    """
    Runs all configuration validations and returns collected errors.
    This is intended to be called after check_config() and basic setup,
    but before pipeline steps begin.

    Parameters:
    config (dict): The loaded configuration dictionary.
    mode (str): "bulk" for Shiba/SnakeShiba, "sc" for scShiba/SnakeScShiba.

    Returns:
    list: A list of all error message strings (empty if everything is valid).
    """
    errors = []

    # 1. Validate that essential files exist
    if config.get('gtf'):
        errors.extend(validate_file_exists(config['gtf'], 'GTF file'))
    if config.get('experiment_table'):
        errors.extend(validate_file_exists(config['experiment_table'], 'Experiment table'))

    # 2. Validate experiment table columns
    if config.get('experiment_table') and os.path.isfile(config['experiment_table']):
        if mode == "bulk":
            errors.extend(
                validate_experiment_table_columns(
                    config['experiment_table'], ["sample", "bam", "group"]
                )
            )
        elif mode == "sc":
            errors.extend(
                validate_experiment_table_columns(
                    config['experiment_table'], ["barcode", "SJ"]
                )
            )

    # 3. Validate reference_group and alternative_group
    only_psi = config.get('only_psi', False)
    only_psi_group = config.get('only_psi_group', False)
    ref_group = config.get('reference_group')
    alt_group = config.get('alternative_group')

    if not (only_psi and only_psi_group):
        # Check that reference_group and alternative_group are specified in config
        if not ref_group:
            errors.append('reference_group is not specified in configuration file')
        if not alt_group:
            errors.append('alternative_group is not specified in configuration file')

        # Validate group names against experiment table
        if ref_group and alt_group and config.get('experiment_table') and os.path.isfile(config['experiment_table']):
            if mode == "bulk":
                errors.extend(
                    validate_groups_bulk(config['experiment_table'], ref_group, alt_group)
                )
            elif mode == "sc":
                errors.extend(
                    validate_groups_sc(config['experiment_table'], ref_group, alt_group)
                )

    # 4. Validate config parameter types and ranges
    errors.extend(validate_config_types(config, mode=mode))

    return errors
