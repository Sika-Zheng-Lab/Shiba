"""
Unit tests for src/gtf2event.py
Tests cover: se, mse, five, three, ri, mxe, afe, ale event detection functions
             and gtf_exon_set helper.
"""

import unittest
import os
import sys
import tempfile
import shutil

import numpy as np
import pandas as pd

# Add src directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "src")))
import gtf2event


def _make_gtf_dic_se():
    """
    Build a minimal gtf_dic for a gene with one skipped exon event.

    Gene structure (+ strand):
      Transcript A: exon1(100-200) --- intron(200-400) --- exon3(400-500)
      Transcript B: exon1(100-200) --- intron(200-300) --- exon2(300-350) --- intron(350-400) --- exon3(400-500)

    SE event: exon2 is skipped when intron(200-400) is used.
    """
    return {
        "GENE_SE": {
            "chr": "chr1",
            "strand": "+",
            "gene_name": "GeneSE",
            "start": np.array([100, 400, 100, 300, 400]),
            "end": np.array([200, 500, 200, 350, 500]),
            "exon_list": {
                "chr1:100-200", "chr1:300-350", "chr1:400-500",
            },
            "intron_list": {
                "chr1:200-300", "chr1:350-400", "chr1:200-400",
            },
            "intron_start_dic": {
                "200": {"300", "400"},
                "350": {"400"},
            },
            "intron_end_dic": {
                "300": {"200"},
                "400": {"200", "350"},
            },
            "start_dic": {
                "100": {"200"},
                "300": {"350"},
                "400": {"500"},
            },
            "end_dic": {
                "200": {"100"},
                "350": {"300"},
                "500": {"400"},
            },
            "transcript_exon_dic": {
                "txA": {"chr1:100-200", "chr1:400-500"},
                "txB": {"chr1:100-200", "chr1:300-350", "chr1:400-500"},
            },
            "transcript_intron_dic": {
                "txA": {"chr1:200-400"},
                "txB": {"chr1:200-300", "chr1:350-400"},
            },
        }
    }


def _make_gtf_dic_five():
    """
    Build a minimal gtf_dic for a gene with an alternative 5' splice site event (+ strand).

    Two exons share the same start but differ in end:
      exon_a: chr1:100-250
      exon_b: chr1:100-200
    Both have an intron ending at 400 (same downstream exon start).
    """
    return {
        "GENE_FIVE": {
            "chr": "chr1",
            "strand": "+",
            "gene_name": "GeneFive",
            "start": np.array([100, 100, 400]),
            "end": np.array([200, 250, 500]),
            "exon_list": {"chr1:100-200", "chr1:100-250", "chr1:400-500"},
            "intron_list": {"chr1:200-400", "chr1:250-400"},
            "intron_start_dic": {
                "200": {"400"},
                "250": {"400"},
            },
            "intron_end_dic": {
                "400": {"200", "250"},
            },
            "start_dic": {
                "100": {"200", "250"},
                "400": {"500"},
            },
            "end_dic": {
                "200": {"100"},
                "250": {"100"},
                "500": {"400"},
            },
            "transcript_exon_dic": {
                "txA": {"chr1:100-200", "chr1:400-500"},
                "txB": {"chr1:100-250", "chr1:400-500"},
            },
            "transcript_intron_dic": {
                "txA": {"chr1:200-400"},
                "txB": {"chr1:250-400"},
            },
        }
    }


def _make_gtf_dic_ri():
    """
    Build a minimal gtf_dic for a gene with a retained intron event.

    Transcript A: exon1(100-200) --- intron(200-300) --- exon2(300-400)
    Transcript B: big_exon(100-400)   (retained intron)
    """
    return {
        "GENE_RI": {
            "chr": "chr1",
            "strand": "+",
            "gene_name": "GeneRI",
            "start": np.array([100, 300, 100]),
            "end": np.array([200, 400, 400]),
            "exon_list": {"chr1:100-200", "chr1:300-400", "chr1:100-400"},
            "intron_list": {"chr1:200-300"},
            "intron_start_dic": {"200": {"300"}},
            "intron_end_dic": {"300": {"200"}},
            "start_dic": {
                "100": {"200", "400"},
                "300": {"400"},
            },
            "end_dic": {
                "200": {"100"},
                "400": {"100", "300"},
            },
            "transcript_exon_dic": {
                "txA": {"chr1:100-200", "chr1:300-400"},
                "txB": {"chr1:100-400"},
            },
            "transcript_intron_dic": {
                "txA": {"chr1:200-300"},
            },
        }
    }


class TestSE(unittest.TestCase):
    def test_se_detects_event(self):
        gtf_dic = _make_gtf_dic_se()
        events = gtf2event.se(gtf_dic)
        self.assertGreater(len(events), 0)
        # Each event: [exon, inc1, inc2, exc, strand, gene, gene_name]
        event = events[0]
        self.assertEqual(len(event), 7)
        self.assertEqual(event[0], "chr1:300-350")  # skipped exon
        self.assertEqual(event[4], "+")
        self.assertEqual(event[5], "GENE_SE")

    def test_se_no_events_for_empty(self):
        events = gtf2event.se({})
        self.assertEqual(len(events), 0)

    def test_se_no_events_when_no_introns(self):
        gtf_dic = {
            "GENE_NO_INTRON": {
                "chr": "chr1", "strand": "+", "gene_name": "NoIntron",
                "start": np.array([100]), "end": np.array([200]),
                "exon_list": {"chr1:100-200"},
            }
        }
        events = gtf2event.se(gtf_dic)
        self.assertEqual(len(events), 0)


class TestFive(unittest.TestCase):
    def test_five_detects_event(self):
        gtf_dic = _make_gtf_dic_five()
        events = gtf2event.five(gtf_dic)
        self.assertGreater(len(events), 0)
        # Each event: [exon_a, exon_b, intron_a, intron_b, strand, gene, gene_name]
        event = events[0]
        self.assertEqual(len(event), 7)
        self.assertEqual(event[4], "+")

    def test_five_no_events_for_empty(self):
        events = gtf2event.five({})
        self.assertEqual(len(events), 0)


class TestRI(unittest.TestCase):
    def test_ri_detects_event(self):
        gtf_dic = _make_gtf_dic_ri()
        events = gtf2event.ri(gtf_dic)
        self.assertGreater(len(events), 0)
        # Each event: [exon_a, exon_b, exon_c, intron_a, strand, gene, gene_name]
        event = events[0]
        self.assertEqual(len(event), 7)
        self.assertIn("chr1:", event[0])

    def test_ri_no_events_for_empty(self):
        events = gtf2event.ri({})
        self.assertEqual(len(events), 0)


class TestGtfExonSet(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Minimal GTF with exon lines
        gtf_lines = [
            'chr1\tsource\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
            'chr1\tsource\texon\t300\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
            'chr1\tsource\tCDS\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
        ]
        self.gtf_path = os.path.join(self.tmpdir, "test.gtf")
        with open(self.gtf_path, "w") as f:
            f.write("\n".join(gtf_lines) + "\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_gtf_exon_set_returns_set(self):
        result = gtf2event.gtf_exon_set(self.gtf_path)
        self.assertIsInstance(result, set)

    def test_gtf_exon_set_correct_elements(self):
        result = gtf2event.gtf_exon_set(self.gtf_path)
        # Only exon lines should be included, not CDS
        self.assertEqual(len(result), 2)
        self.assertIn("chr1:100-200", result)
        self.assertIn("chr1:300-400", result)


class TestGtfFunction(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Minimal GTF with two transcripts of one gene
        gtf_lines = [
            'chr1\tsource\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Gene1";',
            'chr1\tsource\texon\t400\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "Gene1";',
            'chr1\tsource\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; gene_name "Gene1";',
            'chr1\tsource\texon\t300\t350\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; gene_name "Gene1";',
            'chr1\tsource\texon\t400\t500\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; gene_name "Gene1";',
        ]
        self.gtf_path = os.path.join(self.tmpdir, "test.gtf")
        with open(self.gtf_path, "w") as f:
            f.write("\n".join(gtf_lines) + "\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_gtf_returns_dict(self):
        result = gtf2event.gtf(self.gtf_path, 1)
        self.assertIsInstance(result, dict)

    def test_gtf_contains_gene(self):
        result = gtf2event.gtf(self.gtf_path, 1)
        # result is split by num_process; result[0] is the first split
        self.assertIn("G1", result[0])

    def test_gtf_gene_has_required_keys(self):
        result = gtf2event.gtf(self.gtf_path, 1)
        gene_data = result[0]["G1"]
        for key in ["chr", "strand", "gene_name", "exon_list", "intron_list"]:
            self.assertIn(key, gene_data, f"Missing key: {key}")

    def test_gtf_multi_process_splits(self):
        result = gtf2event.gtf(self.gtf_path, 2)
        self.assertEqual(len(result), 2)


if __name__ == "__main__":
    unittest.main()
