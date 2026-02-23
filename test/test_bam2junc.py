"""
Unit tests for src/bam2junc.py
Tests cover: prepare_output_dir, create_saf_file, merge_junction_files
(process_samples is skipped since it requires external tools: regtools, featureCounts)
"""

import unittest
import os
import sys
import tempfile
import shutil

import pandas as pd

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
import bam2junc

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


class TestPrepareOutputDir(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_creates_directories(self):
        output_path = os.path.join(self.tmpdir, "output", "junctions.bed")
        output_dir, logs_dir, tmp_dir = bam2junc.prepare_output_dir(output_path)
        self.assertTrue(os.path.isdir(logs_dir))
        self.assertTrue(os.path.isdir(tmp_dir))
        self.assertEqual(os.path.basename(logs_dir), "logs")
        self.assertEqual(os.path.basename(tmp_dir), "tmp")


class TestCreateSafFile(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a minimal RI event file
        self.ri_event_path = os.path.join(self.tmpdir, "EVENT_RI.txt")
        with open(self.ri_event_path, "w") as f:
            f.write("event_id\tpos_id\texon_a\texon_b\texon_c\tintron_a\tstrand\tgene_id\tgene_name\tlabel\n")
            f.write("RI_1\tp1\tea\teb\tec\tchr1:15100-15200\t+\tG1\tG1\tannotated\n")
            f.write("RI_2\tp2\tea\teb\tec\tchr2:4100-4200\t-\tG2\tG2\tunannotated\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_create_saf_file(self):
        saf_path = bam2junc.create_saf_file(self.ri_event_path, self.tmpdir)
        self.assertTrue(os.path.isfile(saf_path))
        df = pd.read_csv(saf_path, sep="\t")
        self.assertIn("GeneID", df.columns)
        self.assertIn("Chr", df.columns)
        self.assertIn("Start", df.columns)
        self.assertIn("End", df.columns)
        self.assertIn("Strand", df.columns)
        # 2 RI events → 4 SAF entries (start+end for each), deduplicated
        self.assertEqual(len(df), 4)

    def test_saf_file_content(self):
        saf_path = bam2junc.create_saf_file(self.ri_event_path, self.tmpdir)
        df = pd.read_csv(saf_path, sep="\t")
        # chr1:15100-15200 → start junc chr1:15100-15101, end junc chr1:15199-15200
        gene_ids = set(df["GeneID"].values)
        self.assertIn("chr1:15100-15101", gene_ids)
        self.assertIn("chr1:15199-15200", gene_ids)


class TestMergeJunctionFiles(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create mock exon-exon junction files (regtools output format)
        # regtools format: chr start end name count strand thickStart thickEnd color blockCount blockSizes blockStarts
        self.junc1 = os.path.join(self.tmpdir, "s1_exon-exon.junc")
        with open(self.junc1, "w") as f:
            f.write("chr1\t100\t300\tjunc1\t50\t+\t100\t300\t0\t2\t20,25\t0,175\n")
            f.write("chr1\t400\t600\tjunc2\t30\t+\t400\t600\t0\t2\t15,20\t0,180\n")

        self.junc2 = os.path.join(self.tmpdir, "s2_exon-exon.junc")
        with open(self.junc2, "w") as f:
            f.write("chr1\t100\t300\tjunc1\t45\t+\t100\t300\t0\t2\t20,25\t0,175\n")
            f.write("chr1\t400\t600\tjunc2\t35\t+\t400\t600\t0\t2\t15,20\t0,180\n")

        # Create mock exon-intron junction files (featureCounts output format)
        self.ei_junc1 = os.path.join(self.tmpdir, "s1_exon-intron.junc")
        with open(self.ei_junc1, "w") as f:
            f.write("# Program:featureCounts\n")
            f.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\ts1.bam\n")
            f.write("chr1:500-501\tchr1\t500\t501\t+\t1\t20\n")

        self.ei_junc2 = os.path.join(self.tmpdir, "s2_exon-intron.junc")
        with open(self.ei_junc2, "w") as f:
            f.write("# Program:featureCounts\n")
            f.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\ts2.bam\n")
            f.write("chr1:500-501\tchr1\t500\t501\t+\t1\t25\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_merge_junction_files(self):
        junc_files = [
            (self.junc1, "exon-exon"),
            (self.junc2, "exon-exon"),
            (self.ei_junc1, "exon-intron"),
            (self.ei_junc2, "exon-intron"),
        ]
        output_path = os.path.join(self.tmpdir, "merged_junctions.bed")
        bam2junc.merge_junction_files(junc_files, output_path)
        self.assertTrue(os.path.isfile(output_path))
        result = pd.read_csv(output_path, sep="\t")
        self.assertIn("chr", result.columns)
        self.assertIn("ID", result.columns)
        self.assertGreater(len(result), 0)


if __name__ == "__main__":
    unittest.main()
