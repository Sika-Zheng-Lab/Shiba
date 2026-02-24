"""
Unit tests for src/pca.py
Tests cover: load_tpm_table, logit_conversion, load_psi_table, mtx2pca
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
import pca

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


class TestLogitConversion(unittest.TestCase):
    def test_logit_basic(self):
        """logit of 0.5 should be 0."""
        result = pca.logit_conversion(0.5)
        self.assertAlmostEqual(result, 0.0, places=5)

    def test_logit_high_value(self):
        """logit of value near 1 should be large positive."""
        result = pca.logit_conversion(0.99)
        self.assertGreater(result, 0)

    def test_logit_low_value(self):
        """logit of value near 0 should be large negative."""
        result = pca.logit_conversion(0.01)
        self.assertLess(result, 0)

    def test_logit_clamps_zero_and_one(self):
        """Values at 0 and 1 should be clamped to epsilon boundaries."""
        result_zero = pca.logit_conversion(0.0)
        result_one = pca.logit_conversion(1.0)
        self.assertTrue(np.isfinite(result_zero))
        self.assertTrue(np.isfinite(result_one))


class TestLoadTpmTable(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a small TPM table
        df = pd.DataFrame({
            "Gene": ["G1", "G2", "G3"],
            "s1": [10.0, 0.0, 5.0],
            "s2": [8.0, 0.0, 3.0],
            "s3": [12.0, 0.0, 7.0],
        })
        self.tpm_path = os.path.join(self.tmpdir, "tpm.tsv")
        df.to_csv(self.tpm_path, sep="\t", index=False)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_load_tpm_table_drops_zero_rows(self):
        result = pca.load_tpm_table(self.tpm_path)
        # G2 has all zeros, should be dropped; Gene is the index (index_col=0)
        self.assertEqual(result.shape[0], 2)
        self.assertNotIn("G2", result.index)

    def test_load_tpm_table_preserves_values(self):
        result = pca.load_tpm_table(self.tpm_path)
        # load_tpm_table does NOT log-transform; values remain as-is
        self.assertAlmostEqual(result.loc["G1", "s1"], 10.0, places=2)


class TestLoadPsiTable(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a small PSI table with some NaN
        df = pd.DataFrame({
            "event_id": ["E1", "E2", "E3"],
            "pos_id": ["p1", "p2", "p3"],
            "s1": [0.8, np.nan, 0.5],
            "s2": [0.7, 0.6, np.nan],
            "s3": [0.9, 0.8, 0.4],
            "s4": [0.85, 0.7, 0.45],
        })
        self.psi_path = os.path.join(self.tmpdir, "psi.tsv")
        df.to_csv(self.psi_path, sep="\t", index=False)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_load_psi_table_shape(self):
        result = pca.load_psi_table(self.psi_path)
        # event_id is the index (index_col=0), pos_id is dropped
        self.assertEqual(result.index.name, "event_id")
        self.assertNotIn("pos_id", result.columns)
        self.assertGreater(result.shape[0], 0)

    def test_load_psi_table_no_nan_after_imputation(self):
        result = pca.load_psi_table(self.psi_path)
        # After KNN imputation, no NaN should remain; all cols are samples
        self.assertFalse(result.isna().any().any())


class TestMtx2Pca(unittest.TestCase):
    def test_mtx2pca_output(self):
        """mtx2pca should return (pca_df, contribution_df)."""
        df = pd.DataFrame({
            "s1": [1.0, 2.0, 3.0, 4.0, 5.0],
            "s2": [1.1, 2.1, 3.1, 4.1, 5.1],
            "s3": [5.0, 4.0, 3.0, 2.0, 1.0],
            "s4": [5.1, 4.1, 3.1, 2.1, 1.1],
        }, index=["G1", "G2", "G3", "G4", "G5"])
        pca_df, contrib_df = pca.mtx2pca(df, genes=3000)
        self.assertEqual(pca_df.shape[0], 4)  # 4 samples
        self.assertGreater(pca_df.shape[1], 0)

    def test_mtx2pca_fewer_samples(self):
        """When fewer samples than features, PCA should still work."""
        df = pd.DataFrame({
            "s1": [1.0, 2.0, 3.0],
            "s2": [3.0, 4.0, 1.0],
        }, index=["G1", "G2", "G3"])
        pca_df, contrib_df = pca.mtx2pca(df, genes=3000)
        self.assertEqual(pca_df.shape[0], 2)


if __name__ == "__main__":
    unittest.main()
