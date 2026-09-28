"""Unit tests for Tosa single-cell output conversion."""

import csv
import gzip
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
import sc2junc


class TestSingleCellTosaAdapter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.ri = os.path.join(self.tmp.name, "EVENT_RI.txt")
        with open(self.ri, "w") as handle:
            handle.write("event_id\tintron_a\nRI_1\tchr1:100-200\nRI_2\tchr2:300-400\n")

    def barcode_table(self, name, rows):
        path = os.path.join(self.tmp.name, name + "_barcodes.tsv")
        with open(path, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["barcode", "group"])
            writer.writerows(rows)
        return path

    def tosa_output(self, name, junctions, boundaries):
        prefix = os.path.join(self.tmp.name, name)
        with gzip.open(prefix + "_junction_barcodes.tsv.gz", "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Feature", "Strand", "Barcode", "Count"])
            writer.writerows(junctions)
        with gzip.open(prefix + "_boundary_barcodes_detail.tsv.gz", "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Boundary", "Type", "Strand", "Barcode", "Count"])
            writer.writerows(boundaries)
        return prefix

    def test_grouping_across_libraries_strands_and_ri_zeros(self):
        first = self.tosa_output("s1", [["chr1:101-199", "+", "A", 2],
                                         ["chr1:101-199", "-", "B", 3],
                                         ["chr1:101-199", "+", "X", 99]],
                                 [["chr1:99-101", "5p", "+", "A", 1],
                                  ["chr1:198-200", "3p", "+", "B", 2],
                                  ["chr1:499-501", "5p", "+", "A", 8]])
        second = self.tosa_output("s2", [["chr1:101-199", "+", "A", 4]],
                                  [["chr1:99-101", "5p", "+", "A", 5]])
        samples = {
            "s1": {"barcode": self.barcode_table("s1", [["A", "Ref"], ["B", "Alt"]])},
            "s2": {"barcode": self.barcode_table("s2", [["A", "Alt"]])},
        }
        output = os.path.join(self.tmp.name, "junctions.bed")
        sc2junc.merge_tosa_samples(samples, {"s1": first, "s2": second}, self.ri, output)
        with open(output) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, ["chr", "start", "end", "ID", "Alt", "Ref"])
            rows = {row["ID"]: row for row in reader}
        self.assertEqual(len(rows), 5)
        self.assertEqual((rows["chr1:100-200"]["Ref"], rows["chr1:100-200"]["Alt"]), ("2", "7"))
        self.assertEqual((rows["chr1:100-101"]["Ref"], rows["chr1:100-101"]["Alt"]), ("1", "5"))
        self.assertEqual(rows["chr1:199-200"]["Alt"], "2")
        self.assertEqual(rows["chr2:300-301"]["Ref"], "0")
        self.assertNotIn("chr1:500-501", rows)

    def test_boundary_anchor_two_and_strand_option(self):
        prefix = self.tosa_output("s", [], [["chr1:98-102", "5p", "+", "A", 3],
                                           ["chr1:197-201", "3p", "+", "A", 4]])
        counts = {"Ref": sc2junc.Counter()}
        sc2junc.add_tosa_counts(prefix, {"A": "Ref"}, 2, sc2junc.ri_boundary_ids(self.ri), counts)
        self.assertEqual(counts["Ref"]["chr1:100-101"], 3)
        self.assertEqual(counts["Ref"]["chr1:199-200"], 4)
        args = ("x.bam", "x.gtf", "out", "barcodes.tsv", 2, 6, 1, 70, 500000)
        self.assertNotIn("-s", sc2junc.tosa_command(*args, "unstranded"))
        self.assertIn("RF", sc2junc.tosa_command(*args, "1"))

    def test_conflicting_barcode_groups_rejected(self):
        path = self.barcode_table("s", [["A", "Ref"], ["A", "Alt"]])
        with self.assertRaisesRegex(ValueError, "conflicting groups"):
            sc2junc.read_barcode_groups(path)


if __name__ == "__main__":
    unittest.main()
