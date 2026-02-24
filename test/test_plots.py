"""
Unit tests for src/plots.py
Tests cover: load_experiment_table, calculate_event_count, get_shiba_version,
             create_empty_plot_content, load_splicing_summary_table
"""

import unittest
import os
import sys
import tempfile
import shutil

import pandas as pd

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
import plots


class TestLoadExperimentTable(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.exp_path = os.path.join(self.tmpdir, "experiment.tsv")
        with open(self.exp_path, "w") as f:
            f.write("sample\tbam\tgroup\n")
            f.write("s1\ts1.bam\tctrl\n")
            f.write("s2\ts2.bam\tctrl\n")
            f.write("s3\ts3.bam\ttreat\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_returns_dataframe(self):
        result = plots.load_experiment_table(self.exp_path)
        self.assertIsInstance(result, pd.DataFrame)

    def test_has_group_order_column(self):
        result = plots.load_experiment_table(self.exp_path)
        self.assertIn("group_order", result.columns)

    def test_group_order_values(self):
        result = plots.load_experiment_table(self.exp_path)
        # ctrl appears first -> 0, treat -> 1
        ctrl_order = result.loc[result["group"] == "ctrl", "group_order"].iloc[0]
        treat_order = result.loc[result["group"] == "treat", "group_order"].iloc[0]
        self.assertEqual(ctrl_order, 0)
        self.assertEqual(treat_order, 1)

    def test_correct_sample_count(self):
        result = plots.load_experiment_table(self.exp_path)
        self.assertEqual(len(result), 3)


class TestCalculateEventCount(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create input_dir/splicing/PSI_SE.txt
        splicing_dir = os.path.join(self.tmpdir, "splicing")
        os.makedirs(splicing_dir)
        df = pd.DataFrame({
            "event_id": ["E1", "E2", "E3"],
            "dPSI": [0.1, -0.2, 0.05],
        })
        df.to_csv(os.path.join(splicing_dir, "PSI_SE.txt"), sep="\t", index=False)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_count_events(self):
        count = plots.calculate_event_count(self.tmpdir, "SE")
        self.assertEqual(count, 3)

    def test_missing_file_returns_zero(self):
        count = plots.calculate_event_count(self.tmpdir, "NONEXISTENT")
        self.assertEqual(count, 0)


class TestGetShibaVersion(unittest.TestCase):
    def test_returns_string(self):
        version = plots.get_shiba_version()
        self.assertIsInstance(version, str)


class TestCreateEmptyPlotContent(unittest.TestCase):
    def test_returns_html_with_message(self):
        result = plots.create_empty_plot_content("No data")
        self.assertIn("No data", result)
        self.assertIn("<div", result)

    def test_escapes_html_chars(self):
        result = plots.create_empty_plot_content("<script>alert('xss')</script>")
        self.assertNotIn("<script>alert", result)


class TestLoadSplicingSummaryTable(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        splicing_dir = os.path.join(self.tmpdir, "splicing")
        os.makedirs(splicing_dir)
        df = pd.DataFrame({
            "AS": ["SE", "SE", "FIVE", "FIVE"],
            "Direction": ["up", "down", "up", "down"],
            "Number": [10, 5, 3, 2],
        })
        df.to_csv(os.path.join(splicing_dir, "summary.txt"), sep="\t", index=False)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_returns_dataframe(self):
        result = plots.load_splicing_summary_table(self.tmpdir)
        self.assertIsInstance(result, pd.DataFrame)

    def test_correct_shape(self):
        result = plots.load_splicing_summary_table(self.tmpdir)
        self.assertEqual(result.shape[0], 4)
        self.assertIn("AS", result.columns)


if __name__ == "__main__":
    unittest.main()
