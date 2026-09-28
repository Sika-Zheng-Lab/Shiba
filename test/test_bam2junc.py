"""Tests for Tosa output conversion into Shiba's junction matrix."""

import csv
import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
import bam2junc


class TestTosaAdapter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.ri = os.path.join(self.tmp.name, "EVENT_RI.txt")
        with open(self.ri, "w") as handle:
            handle.write("event_id\tintron_a\nRI_1\tchr1:100-200\nRI_2\tchr2:300-400\n")

    def write_tosa(self, sample, junction_rows, boundary_rows):
        prefix = os.path.join(self.tmp.name, sample)
        with gzip.open(prefix + "_junction.tsv.gz", "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Junction", "Strand", "Count"])
            writer.writerows(junction_rows)
        with gzip.open(prefix + "_boundary.tsv.gz", "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Boundary", "Type", "Strand", "Count"])
            writer.writerows(boundary_rows)
        return prefix

    def test_coordinates_strands_ri_zeros_and_sample_names(self):
        second = self.write_tosa("sample_2", [["1:101-199", "+", 3]],
                                 [["1:99-101", "5p", "+", 2],
                                  ["1:198-200", "3p", "+", 4]])
        first = self.write_tosa("sample_1", [["1:101-199", "+", 5],
                                              ["1:101-199", "-", 7]],
                                [["1:99-101", "5p", "+", 6],
                                 ["1:99-101", "5p", "-", 1],
                                 ["1:198-200", "3p", "+", 8],
                                 ["1:500-502", "5p", "+", 9]])
        output = os.path.join(self.tmp.name, "junctions.bed")
        bam2junc.merge_tosa_samples({"sample_2": second, "sample_1": first}, self.ri, output)
        with open(output) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, ["chr", "start", "end", "ID", "sample_1", "sample_2"])
            rows = {row["ID"]: row for row in reader}
        self.assertEqual(len(rows), 5)
        self.assertEqual((rows["chr1:100-200"]["sample_1"], rows["chr1:100-200"]["sample_2"]), ("12", "3"))
        self.assertEqual((rows["chr1:100-101"]["sample_1"], rows["chr1:100-101"]["sample_2"]), ("7", "2"))
        self.assertEqual(rows["chr1:199-200"]["sample_1"], "8")
        self.assertEqual(rows["chr2:300-301"]["sample_1"], "0")
        self.assertNotIn("chr1:501-502", rows)

    def test_boundary_anchor_two(self):
        prefix = self.write_tosa("s", [], [["chr1:98-102", "5p", "+", 3],
                                                  ["chr1:197-201", "3p", "+", 4]])
        counts = bam2junc.read_tosa_counts(prefix, 2, bam2junc.ri_boundary_ids(self.ri))
        self.assertEqual(counts["chr1:100-101"], 3)
        self.assertEqual(counts["chr1:199-200"], 4)

    def test_strand_translation(self):
        args = ("x.bam", "x.gtf", "out", 2, 6, 1, 70, 500000)
        self.assertIn("RF", bam2junc.tosa_command(*args, "1"))
        self.assertNotIn("-s", bam2junc.tosa_command(*args, "unstranded"))
        with self.assertRaises(ValueError):
            bam2junc.tosa_command(*args, "invalid")


if __name__ == "__main__":
    unittest.main()
