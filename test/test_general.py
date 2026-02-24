import unittest
from unittest.mock import patch, mock_open, MagicMock
import subprocess
import os
import sys
import yaml
import json
import datetime
import tempfile
import shutil

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
from lib.general import load_config, check_config, execute_command, generate_report, check_samplesize, check_groupsize

class TestGeneralModule(unittest.TestCase):
    def test_load_config_success(self):
        mock_yaml_content = """
        key1: value1
        key2: value2
        """
        with patch("builtins.open", mock_open(read_data=mock_yaml_content)):
            with patch("yaml.safe_load", return_value=yaml.safe_load(mock_yaml_content)) as mock_safe_load:
                config = load_config("config.yaml")
                self.assertEqual(config, {'key1': 'value1', 'key2': 'value2'})
                mock_safe_load.assert_called_once()

    def test_load_config_file_not_found(self):
        with patch("builtins.open", side_effect=FileNotFoundError):
            with self.assertRaises(FileNotFoundError):
                load_config("nonexistent_config.yaml")

    def test_load_config_yaml_error(self):
        invalid_yaml_content = "{invalid: yaml,"
        with patch("builtins.open", mock_open(read_data=invalid_yaml_content)):
            with patch("yaml.safe_load", side_effect=yaml.YAMLError):
                with self.assertRaises(ValueError):
                    load_config("invalid_config.yaml")

    def test_check_config_all_keys_present(self):
        config = {'key1': 'value1', 'key2': 'value2'}
        keys = ['key1', 'key2']
        missing_keys = check_config(config, keys)
        self.assertEqual(missing_keys, [])

    def test_check_config_missing_keys(self):
        config = {'key1': 'value1'}
        keys = ['key1', 'key2']
        missing_keys = check_config(config, keys)
        self.assertEqual(missing_keys, ['key2'])

    def test_check_config_empty_values(self):
        config = {'key1': 'value1', 'key2': None}
        keys = ['key1', 'key2']
        missing_keys = check_config(config, keys)
        self.assertEqual(missing_keys, ['key2'])

    @patch("subprocess.run")
    def test_execute_command_success(self, mock_subprocess):
        mock_subprocess.return_value = MagicMock(returncode=0)
        returncode = execute_command(["echo", "Hello"])
        self.assertEqual(returncode, 0)
        mock_subprocess.assert_called_once_with(["echo", "Hello"], shell=False)

    @patch("subprocess.run")
    def test_execute_command_failure(self, mock_subprocess):
        mock_subprocess.return_value = MagicMock(returncode=1)
        returncode = execute_command(["false"])
        self.assertEqual(returncode, 1)
        mock_subprocess.assert_called_once_with(["false"], shell=False)

    @patch("subprocess.run")
    def test_execute_command_with_log(self, mock_subprocess):
        mock_subprocess.return_value = MagicMock(returncode=0)
        tmpdir = tempfile.mkdtemp()
        try:
            log_file = os.path.join(tmpdir, "test.log")
            returncode = execute_command(["echo", "test"], log_file)
            self.assertEqual(returncode, 0)
        finally:
            shutil.rmtree(tmpdir)


class TestGenerateReport(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_generate_report_creates_json(self):
        start = datetime.datetime.now()
        report_path = generate_report(
            name="Shiba",
            output_dir=self.tmpdir,
            version="1.0.0",
            command_line="shiba run config.yaml",
            experiment_table="/path/to/exp.tsv",
            start_time=start
        )
        self.assertTrue(os.path.isfile(report_path))
        self.assertTrue(report_path.endswith("report.json"))

    def test_generate_report_content(self):
        start = datetime.datetime(2024, 1, 1, 12, 0, 0)
        report_path = generate_report(
            name="Shiba",
            output_dir=self.tmpdir,
            version="2.0.0",
            command_line="shiba run config.yaml",
            experiment_table="/path/to/exp.tsv",
            start_time=start
        )
        with open(report_path, "r") as f:
            data = json.load(f)
        self.assertEqual(data["tool"]["name"], "Shiba")
        self.assertEqual(data["tool"]["version"], "2.0.0")
        self.assertEqual(data["run"]["command"], "shiba run config.yaml")
        self.assertEqual(data["experiment"]["table_path"], "/path/to/exp.tsv")
        self.assertIn("start_time", data["run"])
        self.assertIn("duration_seconds", data["run"])

    def test_generate_report_no_start_time(self):
        report_path = generate_report(
            name="Shiba",
            output_dir=self.tmpdir,
            version="1.0.0",
            command_line="shiba run config.yaml",
            experiment_table="/path/to/exp.tsv",
            start_time=None
        )
        with open(report_path, "r") as f:
            data = json.load(f)
        self.assertEqual(data["run"]["start_time"], "NA")
        self.assertEqual(data["run"]["duration_seconds"], "NA")


class TestCheckSamplesize(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_check_samplesize_valid(self):
        exp_path = os.path.join(self.tmpdir, "exp.tsv")
        with open(exp_path, "w") as f:
            f.write("sample\tbam\tgroup\n")
            f.write("s1\ts1.bam\tctrl\n")
            f.write("s2\ts2.bam\ttreat\n")
        result = check_samplesize(exp_path)
        self.assertEqual(result, 2)

    def test_check_samplesize_empty(self):
        exp_path = os.path.join(self.tmpdir, "exp.tsv")
        with open(exp_path, "w") as f:
            f.write("sample\tbam\tgroup\n")
        with self.assertRaises(ValueError):
            check_samplesize(exp_path)

    def test_check_samplesize_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            check_samplesize("/nonexistent/path.tsv")


class TestCheckGroupsize(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_check_groupsize_valid(self):
        exp_path = os.path.join(self.tmpdir, "exp.tsv")
        with open(exp_path, "w") as f:
            f.write("sample\tbam\tgroup\n")
            f.write("s1\ts1.bam\tctrl\n")
            f.write("s2\ts2.bam\tctrl\n")
            f.write("s3\ts3.bam\ttreat\n")
        result = check_groupsize(exp_path)
        self.assertEqual(result, 2)

    def test_check_groupsize_single_group(self):
        exp_path = os.path.join(self.tmpdir, "exp.tsv")
        with open(exp_path, "w") as f:
            f.write("sample\tbam\tgroup\n")
            f.write("s1\ts1.bam\tctrl\n")
        result = check_groupsize(exp_path)
        self.assertEqual(result, 1)

    def test_check_groupsize_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            check_groupsize("/nonexistent/path.tsv")

if __name__ == "__main__":
    unittest.main()
