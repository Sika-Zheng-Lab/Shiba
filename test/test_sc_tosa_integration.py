"""Run scShiba on generated tagged BAMs when Tosa and pysam are installed."""

import csv
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


@unittest.skipUnless(shutil.which("tosa") and importlib.util.find_spec("pysam"),
                     "Tosa and pysam are required for the BAM integration test")
class TestSyntheticSingleCellPipeline(unittest.TestCase):
    def test_group_counts_and_ri_psi(self):
        sys.path.insert(0, str(ROOT / "test"))
        from make_sc_tosa_fixture import build_fixture

        with tempfile.TemporaryDirectory() as tmp:
            directory = build_fixture(tmp)
            result = subprocess.run(
                [sys.executable, str(ROOT / "scshiba.py"), str(directory / "config.yaml"), "-p", "2"],
                cwd=ROOT, env=os.environ.copy(), capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr[-3000:])
            output = directory / "scshiba"
            with open(output / "junctions/junctions.bed") as handle:
                counts = {row["ID"]: row for row in csv.DictReader(handle, delimiter="\t")}
            with open(directory / "expected_counts.json") as handle:
                expected = json.load(handle)
            for junction, group_counts in expected.items():
                for group, value in group_counts.items():
                    self.assertEqual(int(counts[junction][group]), value)
            with open(output / "results/PSI_RI.txt") as handle:
                ri = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(ri), 1)
            self.assertTrue(math.isclose(float(ri[0]["ref_PSI"]), 1 / 3))
            self.assertTrue(math.isclose(float(ri[0]["alt_PSI"]), 2 / 5))
            for event in ("SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"):
                self.assertTrue((output / f"results/PSI_{event}.txt").exists())


if __name__ == "__main__":
    unittest.main()
