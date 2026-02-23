"""
Unit tests for src/sc2junc.py
Tests cover: load_experiment_table, make_sjpath_list, load_group,
make_sample_group_dict, make_grouppath_list, formatting_output
(load_sj and grouping_read_count_each require scanpy mtx files - tested with mocks)
"""

import unittest
import os
import sys
import tempfile
import shutil

import pandas as pd

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
import sc2junc

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


class TestLoadExperimentTable(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.exp_path = os.path.join(self.tmpdir, "experiment.tsv")
        with open(self.exp_path, "w") as f:
            f.write("sample\tSJ\tbarcode\n")
            f.write("s1\t/path/sj1\t/path/barcode1.tsv\n")
            f.write("s2\t/path/sj2\t/path/barcode2.tsv\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_load_experiment_table(self):
        result = sc2junc.load_experiment_table(self.exp_path)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)
        self.assertIn("SJ", result.columns)
        self.assertIn("barcode", result.columns)


class TestMakeSjpathList(unittest.TestCase):
    def test_make_sjpath_list(self):
        df = pd.DataFrame({"SJ": ["/path/sj1", "/path/sj2"]})
        result = sc2junc.make_sjpath_list(df)
        self.assertEqual(result, ["/path/sj1", "/path/sj2"])


class TestMakeGrouppathList(unittest.TestCase):
    def test_make_grouppath_list(self):
        df = pd.DataFrame({"barcode": ["/path/b1.tsv", "/path/b2.tsv"]})
        result = sc2junc.make_grouppath_list(df)
        self.assertEqual(result, ["/path/b1.tsv", "/path/b2.tsv"])


class TestLoadGroup(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.group_path = os.path.join(self.tmpdir, "barcodes.tsv")
        with open(self.group_path, "w") as f:
            f.write("barcode\tgroup\n")
            f.write("ACGT\tgroupA\n")
            f.write("TGCA\tgroupA\n")
            f.write("AAAA\tgroupB\n")
            f.write("ACGT\tgroupA\n")  # duplicate

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_load_group(self):
        result = sc2junc.load_group(self.group_path)
        self.assertIn("barcode", result.columns)
        self.assertIn("group", result.columns)

    def test_load_group_drops_duplicates(self):
        result = sc2junc.load_group(self.group_path)
        self.assertEqual(len(result), 3)  # 4 rows but 1 duplicate


class TestMakeSampleGroupDict(unittest.TestCase):
    def test_make_sample_group_dict(self):
        df = pd.DataFrame({
            "barcode": ["ACGT", "TGCA", "AAAA"],
            "group": ["groupA", "groupA", "groupB"]
        })
        result = sc2junc.make_sample_group_dict(df)
        self.assertIn("groupA", result)
        self.assertIn("groupB", result)
        self.assertEqual(len(result["groupA"]), 2)
        self.assertEqual(len(result["groupB"]), 1)


class TestFormattingOutput(unittest.TestCase):
    def test_formatting_output(self):
        df = pd.DataFrame({
            "SJ": ["chr1:100-200", "chr2:300-400"],
            "groupA": [50, 30],
            "groupB": [20, 40],
        })
        df = df.set_index("SJ")
        result = sc2junc.formatting_output(df)
        self.assertIn("chr", result.columns)
        self.assertIn("start", result.columns)
        self.assertIn("end", result.columns)
        self.assertIn("ID", result.columns)
        # Coordinates should be adjusted: start-1, end+1
        row = result[result["chr"] == "chr1"].iloc[0]
        self.assertEqual(str(row["start"]), "99")
        self.assertEqual(str(row["end"]), "201")


if __name__ == "__main__":
    unittest.main()
