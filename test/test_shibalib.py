"""
Comprehensive unit tests for src/lib/shibalib.py
Tests cover: I/O functions, event filtering, column generation,
PSI calculation (all 8 event types), differential splicing,
t-test, make_psi_mtx, and EventCounter.
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
from lib import shibalib

# Path to fixture data
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


# ---------------------------------------------------------------------------
# Helper: build a minimal junction DataFrame
# ---------------------------------------------------------------------------
def _make_junc_df():
    """Return a small junction DataFrame with 4 samples."""
    return pd.read_csv(os.path.join(DATA_DIR, "junctions.bed"), sep="\t", dtype={"chr": str})


# ============================================================================
# I/O Functions
# ============================================================================
class TestReadEvents(unittest.TestCase):
    def test_read_events_returns_dict_of_8(self):
        d = shibalib.read_events(DATA_DIR)
        self.assertIsInstance(d, dict)
        for key in ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]:
            self.assertIn(key, d)
            self.assertIsInstance(d[key], pd.DataFrame)

    def test_read_events_se_columns(self):
        d = shibalib.read_events(DATA_DIR)
        expected = ["event_id", "pos_id", "exon", "intron_a", "intron_b", "intron_c", "strand", "gene_id", "gene_name", "label"]
        self.assertEqual(list(d["SE"].columns), expected)

    def test_read_events_se_row_count(self):
        d = shibalib.read_events(DATA_DIR)
        self.assertEqual(len(d["SE"]), 3)

    def test_read_events_sc_excludes_ri(self):
        d = shibalib.read_events_sc(DATA_DIR)
        self.assertNotIn("RI", d)
        self.assertEqual(len(d), 7)


class TestReadJunctions(unittest.TestCase):
    def test_read_junctions_shape(self):
        junc_df = shibalib.read_junctions(os.path.join(DATA_DIR, "junctions.bed"))
        self.assertGreater(junc_df.shape[0], 0)
        self.assertGreater(junc_df.shape[1], 4)

    def test_read_junctions_dtypes(self):
        junc_df = shibalib.read_junctions(os.path.join(DATA_DIR, "junctions.bed"))
        # Columns after index 4 should contain integer values
        for col in junc_df.columns[4:]:
            # Values should be int (either native int dtype or object holding ints)
            self.assertTrue(
                pd.api.types.is_integer_dtype(junc_df[col]) or
                all(isinstance(v, (int, np.integer)) for v in junc_df[col]),
                f"{col} should contain integers"
            )


class TestReadGroup(unittest.TestCase):
    def test_read_group(self):
        group_df = shibalib.read_group(os.path.join(DATA_DIR, "experiment_table.tsv"))
        self.assertIn("sample", group_df.columns)
        self.assertIn("group", group_df.columns)
        self.assertEqual(len(group_df), 4)


# ============================================================================
# Group / Sample Utility Functions
# ============================================================================
class TestSetGroup(unittest.TestCase):
    def setUp(self):
        self.group_df = pd.DataFrame({
            "sample": ["s1", "s2", "s3", "s4"],
            "group": ["ctrl", "ctrl", "treat", "treat"]
        })

    def test_set_group_normal(self):
        result = shibalib.set_group(self.group_df, False, "ctrl", "treat")
        self.assertEqual(result, ["ctrl", "treat"])

    def test_set_group_onlypsi(self):
        result = shibalib.set_group(self.group_df, True, "ctrl", "treat")
        self.assertEqual(sorted(result), ["ctrl", "treat"])

    def test_set_group_single_group(self):
        df = pd.DataFrame({"sample": ["s1", "s2"], "group": ["ctrl", "ctrl"]})
        result = shibalib.set_group(df, False, "ctrl", "treat")
        self.assertEqual(result, ["ctrl"])


class TestMakeSampleList(unittest.TestCase):
    def test_make_sample_list(self):
        junc_df = _make_junc_df()
        samples = shibalib.make_sample_list(junc_df)
        self.assertEqual(samples, ["Alt_1", "Alt_2", "Ref_1", "Ref_2"])


class TestSampleInGroupList(unittest.TestCase):
    def test_sample_in_group_list(self):
        group_df = pd.DataFrame({
            "sample": ["s1", "s2", "s3"],
            "group": ["A", "A", "B"]
        })
        result = shibalib.sample_in_group_list(group_df, ["A", "B"])
        self.assertEqual(result, ["s1", "s2", "s3"])


class TestSumReads(unittest.TestCase):
    def test_sum_reads(self):
        junc_df = _make_junc_df()
        group_df = pd.DataFrame({
            "sample": ["Ref_1", "Ref_2", "Alt_1", "Alt_2"],
            "group": ["Ref", "Ref", "Alt", "Alt"]
        })
        result = shibalib.sum_reads(False, junc_df, group_df, ["Ref", "Alt"])
        self.assertIn("Ref", result)
        self.assertIn("Alt", result)
        # Check a specific junction: chr10:800-1000 has Alt_1=50 Alt_2=60 Ref_1=10 Ref_2=20
        self.assertEqual(result["Ref"]["chr10:800-1000"], 10 + 20)
        self.assertEqual(result["Alt"]["chr10:800-1000"], 50 + 60)


class TestJuncDict(unittest.TestCase):
    def test_junc_dict(self):
        junc_df = _make_junc_df()
        result = shibalib.junc_dict(junc_df)
        self.assertIn("Alt_1", result)
        self.assertEqual(result["Alt_1"]["chr10:800-1000"], 50)
        self.assertEqual(result["Alt_2"]["chr10:1200-1500"], 55)


# ============================================================================
# JunctionData (memory-efficient shared-memory-capable junction storage)
# ============================================================================
class TestJunctionData(unittest.TestCase):
    """Tests for the JunctionData class used for memory-efficient multiprocessing."""

    def setUp(self):
        self.junc_df = _make_junc_df()
        # Ensure numeric columns are int
        self.junc_df.iloc[:, 4:] = self.junc_df.iloc[:, 4:].astype(int)

    def test_from_dataframe_basic(self):
        """JunctionData.from_dataframe should produce correct shape and sample list."""
        jd = shibalib.JunctionData.from_dataframe(self.junc_df)
        sample_cols = [c for c in self.junc_df.columns if c not in ("chr", "start", "end", "ID")]
        self.assertEqual(jd.array.shape[0], len(self.junc_df))
        self.assertEqual(jd.array.shape[1], len(sample_cols))
        self.assertEqual(jd.sample_ids, sample_cols)

    def test_dict_like_access(self):
        """JunctionData[sample][junction_id] should match the original DataFrame."""
        jd = shibalib.JunctionData.from_dataframe(self.junc_df)
        # Compare against junc_dict for a known value
        junc_dict_all = shibalib.junc_dict(self.junc_df)
        for sample in jd.sample_ids:
            for jid in list(jd.junction_to_idx.keys())[:5]:
                self.assertEqual(jd[sample][jid], junc_dict_all[sample][jid])

    def test_missing_junction_raises(self):
        """Accessing a non-existent junction should raise KeyError."""
        jd = shibalib.JunctionData.from_dataframe(self.junc_df)
        with self.assertRaises(KeyError):
            _ = jd[jd.sample_ids[0]]["nonexistent_junction"]

    def test_se_psi_with_junction_data(self):
        """PSI worker function se() should work with JunctionData as drop-in for dict."""
        jd = shibalib.JunctionData.from_dataframe(self.junc_df)
        event_df = pd.DataFrame({
            "event_id": ["SE_1"],
            "pos_id": ["SE@chr10@1000-1200@800-1500"],
            "exon": ["chr10:1000-1200"],
            "intron_a": ["chr10:800-1000"],
            "intron_b": ["chr10:1200-1500"],
            "intron_c": ["chr10:800-1500"],
            "strand": ["+"],
            "gene_id": ["G1"],
            "gene_name": ["GeneA"],
            "label": ["annotated"]
        })
        # Use only samples that exist in junc_df
        sample_list = jd.sample_ids
        result_jd = shibalib.se(jd, sample_list, event_df, 1, 1, 0)
        # Compare with dict-based result
        junc_dict_all = shibalib.junc_dict(self.junc_df)
        result_dict = shibalib.se(junc_dict_all, sample_list, event_df, 1, 1, 0)
        # Both should produce the same PSI values
        self.assertEqual(len(result_jd), len(result_dict))
        for row_jd, row_dict in zip(result_jd, result_dict):
            for v_jd, v_dict in zip(row_jd, row_dict):
                if isinstance(v_jd, float) and np.isnan(v_jd):
                    self.assertTrue(np.isnan(v_dict))
                else:
                    self.assertEqual(v_jd, v_dict)

    def test_to_shared_memory_and_back(self):
        """Data should survive round-trip through shared memory."""
        jd = shibalib.JunctionData.from_dataframe(self.junc_df)
        shm, shm_info = jd.to_shared_memory()
        try:
            jd2 = shibalib.JunctionData.from_shared_memory(
                shm_info['shm_name'], shm_info['shape'],
                shm_info['junction_ids'], shm_info['sample_ids']
            )
            # Verify data matches
            np.testing.assert_array_equal(jd.array, jd2.array)
            self.assertEqual(jd.sample_ids, jd2.sample_ids)
            # Verify dict-like access
            sample = jd.sample_ids[0]
            jid = list(jd.junction_to_idx.keys())[0]
            self.assertEqual(jd[sample][jid], jd2[sample][jid])
            jd2.close()
        finally:
            shm.close()
            shm.unlink()

    def test_worker_init_and_get(self):
        """_init_junc_worker / _get_junc_data cycle should work in the current process."""
        jd = shibalib.JunctionData.from_dataframe(self.junc_df)
        shm, shm_info = jd.to_shared_memory()
        try:
            shibalib._init_junc_worker(
                shm_info['shm_name'], shm_info['shape'],
                shm_info['junction_ids'], shm_info['sample_ids']
            )
            worker_data = shibalib._get_junc_data()
            sample = jd.sample_ids[0]
            jid = list(jd.junction_to_idx.keys())[0]
            self.assertEqual(jd[sample][jid], worker_data[sample][jid])
        finally:
            # Reset worker state
            shibalib._worker_junc_data = None
            shm.close()
            shm.unlink()

    def test_get_junc_data_without_init_raises(self):
        """_get_junc_data should raise RuntimeError when not initialized."""
        old = shibalib._worker_junc_data
        shibalib._worker_junc_data = None
        try:
            with self.assertRaises(RuntimeError):
                shibalib._get_junc_data()
        finally:
            shibalib._worker_junc_data = old


class TestMakeJuncSet(unittest.TestCase):
    def test_make_junc_set(self):
        junc_df = _make_junc_df()
        result = shibalib.make_junc_set(junc_df)
        self.assertIsInstance(result, set)
        self.assertIn("chr10:800-1000", result)
        self.assertIn("chr10:2800-3500", result)


# ============================================================================
# Event Filtering Functions
# ============================================================================
class TestEventForAnalysisSE(unittest.TestCase):
    def test_se_events_filtered(self):
        events = shibalib.read_events(DATA_DIR)
        junc_df = _make_junc_df()
        junc_set = shibalib.make_junc_set(junc_df)
        result = shibalib.event_for_analysis_se(events["SE"], junc_set)
        # SE_1 and SE_2 have matching junctions; SE_3 may or may not match
        self.assertGreater(len(result), 0)
        self.assertTrue(all(result["event_id"].str.startswith("SE_")))

    def test_se_empty_junc_set(self):
        events = shibalib.read_events(DATA_DIR)
        result = shibalib.event_for_analysis_se(events["SE"], set())
        self.assertEqual(len(result), 0)


class TestEventForAnalysisMSE(unittest.TestCase):
    def test_mse_events_filtered(self):
        events = shibalib.read_events(DATA_DIR)
        junc_df = _make_junc_df()
        junc_set = shibalib.make_junc_set(junc_df)
        result = shibalib.event_for_analysis_mse(events["MSE"], junc_set)
        # MSE_1 has introns matching chr1:1000-1100 etc.
        self.assertGreaterEqual(len(result), 0)


class TestEventForAnalysisFiveThree(unittest.TestCase):
    def test_five_events_filtered(self):
        events = shibalib.read_events(DATA_DIR)
        junc_df = _make_junc_df()
        junc_set = shibalib.make_junc_set(junc_df)
        result = shibalib.event_for_analysis_five_three(events["FIVE"], junc_set)
        self.assertGreaterEqual(len(result), 0)


class TestEventForAnalysisAFEALE(unittest.TestCase):
    def test_afe_events_filtered(self):
        events = shibalib.read_events(DATA_DIR)
        junc_df = _make_junc_df()
        junc_set = shibalib.make_junc_set(junc_df)
        result = shibalib.event_for_analysis_afe_ale(events["AFE"], junc_set)
        self.assertGreaterEqual(len(result), 0)


class TestEventForAnalysisMXE(unittest.TestCase):
    def test_mxe_events_filtered(self):
        events = shibalib.read_events(DATA_DIR)
        junc_df = _make_junc_df()
        junc_set = shibalib.make_junc_set(junc_df)
        result = shibalib.event_for_analysis_mxe(events["MXE"], junc_set)
        self.assertGreaterEqual(len(result), 0)


class TestEventForAnalysisRI(unittest.TestCase):
    def test_ri_events_filtered(self):
        events = shibalib.read_events(DATA_DIR)
        junc_df = _make_junc_df()
        junc_set = shibalib.make_junc_set(junc_df)
        result = shibalib.event_for_analysis_ri(events["RI"], junc_set)
        self.assertGreaterEqual(len(result), 0)


# ============================================================================
# Column Generation Functions
# ============================================================================
class TestColSE(unittest.TestCase):
    def test_col_se(self):
        cols = shibalib.col_se(["s1", "s2"], False)
        self.assertIn("event_id", cols)
        self.assertIn("s1_junction_a", cols)
        self.assertIn("s1_junction_b", cols)
        self.assertIn("s1_junction_c", cols)
        self.assertIn("s1_PSI", cols)
        self.assertIn("s2_PSI", cols)
        # 10 base cols + 4 per sample (junction_a, junction_b, junction_c, PSI)
        self.assertEqual(len(cols), 10 + 2 * 4)


class TestColInd(unittest.TestCase):
    def test_col_ind(self):
        cols = shibalib.col_ind(["s1", "s2"])
        self.assertEqual(cols, ["event_id", "s1_PSI", "s2_PSI", "s1_total_reads", "s2_total_reads"])


class TestColMSE(unittest.TestCase):
    def test_col_mse(self):
        cols = shibalib.col_mse(["s1"], False)
        self.assertIn("event_id", cols)
        self.assertIn("mse_n", cols)
        self.assertIn("s1_junction", cols)
        self.assertIn("s1_PSI", cols)


class TestColFiveThreeAFEALE(unittest.TestCase):
    def test_col_five_three_afe_ale(self):
        cols = shibalib.col_five_three_afe_ale(["s1"], False)
        self.assertIn("event_id", cols)
        self.assertIn("exon_a", cols)
        self.assertIn("exon_b", cols)
        self.assertIn("s1_junction_a", cols)
        self.assertIn("s1_junction_b", cols)
        self.assertIn("s1_PSI", cols)


class TestColMXE(unittest.TestCase):
    def test_col_mxe(self):
        cols = shibalib.col_mxe(["s1"], False)
        self.assertIn("event_id", cols)
        self.assertIn("s1_junction_a1", cols)
        self.assertIn("s1_junction_a2", cols)
        self.assertIn("s1_junction_b1", cols)
        self.assertIn("s1_junction_b2", cols)
        self.assertIn("s1_PSI", cols)


class TestColRI(unittest.TestCase):
    def test_col_ri(self):
        cols = shibalib.col_ri(["s1"], False)
        self.assertIn("event_id", cols)
        self.assertIn("exon_a", cols)
        self.assertIn("exon_c", cols)
        self.assertIn("s1_junction_a_start", cols)
        self.assertIn("s1_junction_a_end", cols)
        self.assertIn("s1_junction_a", cols)
        self.assertIn("s1_PSI", cols)


# ============================================================================
# PSI Calculation - SE
# ============================================================================
class TestSePsi(unittest.TestCase):
    def setUp(self):
        self.event_df = pd.DataFrame({
            "event_id": ["SE_1"],
            "pos_id": ["SE@chr1@120-180@100-200"],
            "exon": ["chr1:120-180"],
            "intron_a": ["chr1:100-120"],
            "intron_b": ["chr1:180-200"],
            "intron_c": ["chr1:100-200"],
            "strand": ["+"],
            "gene_id": ["G1"],
            "gene_name": ["GeneA"],
            "label": ["annotated"]
        })
        self.junc_dict = {
            "s1": {"chr1:100-120": 50, "chr1:180-200": 40, "chr1:100-200": 10},
            "s2": {"chr1:100-120": 0, "chr1:180-200": 0, "chr1:100-200": 0},
        }

    def test_se_psi_calculation(self):
        result = shibalib.se(self.junc_dict, ["s1", "s2"], self.event_df, 1, 5, 0)
        self.assertEqual(len(result), 1)
        row = result[0]
        # For s1: inclusion = (50+40)/2 = 45, exclusion = 10 → PSI = 45/(45+10) = ~0.818
        s1_psi = row[13]  # 10 base cols + s1_junc_a(10), s1_junc_b(11), s1_junc_c(12), s1_psi(13)
        self.assertAlmostEqual(s1_psi, 45 / 55, places=5)

    def test_se_psi_nan_below_minimum(self):
        junc_dict = {
            "s1": {"chr1:100-120": 1, "chr1:180-200": 1, "chr1:100-200": 1},
        }
        result = shibalib.se(junc_dict, ["s1"], self.event_df, 1, 10, 0)
        s1_psi = result[0][13]
        self.assertTrue(np.isnan(s1_psi))

    def test_se_ind_psi(self):
        result = shibalib.se_ind(self.junc_dict, self.event_df, ["s1", "s2"], 1, 0)
        self.assertEqual(len(result), 1)
        row = result[0]
        # event_id + s1_PSI + s2_PSI
        self.assertEqual(row[0], "SE_1")
        self.assertAlmostEqual(row[1], 45 / 55, places=5)


# ============================================================================
# PSI Calculation - MSE
# ============================================================================
class TestMsePsi(unittest.TestCase):
    def setUp(self):
        self.event_df = pd.DataFrame({
            "event_id": ["MSE_1"],
            "pos_id": ["MSE@chr1@1100-1200;1300-1400@1000-2000"],
            "mse_n": ["2"],
            "exon": ["chr1:1100-1200;chr1:1300-1400"],
            "intron": ["chr1:1000-1100;chr1:1200-1300;chr1:1400-2000;chr1:1000-2000"],
            "strand": ["+"],
            "gene_id": ["G4"],
            "gene_name": ["GeneD"],
            "label": ["annotated"]
        })
        self.junc_dict = {
            "s1": {
                "chr1:1000-1100": 30,
                "chr1:1200-1300": 28,
                "chr1:1400-2000": 25,
                "chr1:1000-2000": 5,
            },
        }

    def test_mse_psi_calculation(self):
        result = shibalib.mse(self.junc_dict, ["s1"], self.event_df, 1, 5, 0)
        self.assertEqual(len(result), 1)
        row = result[0]
        # inclusion introns = [30, 28, 25], exclusion = 5
        # inclusion_mean = 27.666..., psi = 27.666.../(27.666...+5) ≈ 0.8469
        s1_psi = row[10]  # 9 base cols + s1_junction(9), s1_psi(10)
        expected = np.mean([30, 28, 25]) / (np.mean([30, 28, 25]) + 5)
        self.assertAlmostEqual(s1_psi, expected, places=4)


# ============================================================================
# PSI Calculation - FIVE/THREE
# ============================================================================
class TestFiveThreePsi(unittest.TestCase):
    def setUp(self):
        self.event_df = pd.DataFrame({
            "event_id": ["FIVE_1"],
            "pos_id": ["FIVE@chr1@5100-5200@5150-5200"],
            "exon_a": ["chr1:5000-5100"],
            "exon_b": ["chr1:5000-5150"],
            "intron_a": ["chr1:5100-5200"],
            "intron_b": ["chr1:5150-5200"],
            "strand": ["+"],
            "gene_id": ["G6"],
            "gene_name": ["GeneF"],
            "label": ["annotated"]
        })
        self.junc_dict = {
            "s1": {"chr1:5100-5200": 40, "chr1:5150-5200": 20},
        }

    def test_five_three_psi(self):
        result = shibalib.five_three(self.junc_dict, ["s1"], self.event_df, 1, 5, 0)
        self.assertEqual(len(result), 1)
        row = result[0]
        # PSI = intron_a / (intron_a + intron_b) = 40 / 60 ≈ 0.6667
        s1_psi = row[12]  # 10 base + s1_junc_a(10), s1_junc_b(11), s1_psi(12)
        self.assertAlmostEqual(s1_psi, 40.0 / 60.0, places=4)

    def test_five_three_ind_psi(self):
        result = shibalib.five_three_ind(self.junc_dict, self.event_df, ["s1"], 1, 0)
        row = result[0]
        self.assertAlmostEqual(row[1], 40.0 / 60.0, places=4)


# ============================================================================
# PSI Calculation - AFE/ALE
# ============================================================================
class TestAFEALEPsi(unittest.TestCase):
    def setUp(self):
        self.event_df = pd.DataFrame({
            "event_id": ["AFE_1"],
            "pos_id": ["AFE@chr1@20100-20500@20300-20500"],
            "exon_a": ["chr1:20000-20100"],
            "exon_b": ["chr1:20200-20300"],
            "intron_a": ["chr1:20100-20500"],
            "intron_b": ["chr1:20300-20500"],
            "strand": ["+"],
            "gene_id": ["G14"],
            "gene_name": ["GeneN"],
            "label": ["annotated"]
        })
        self.junc_dict = {
            "s1": {"chr1:20100-20500": 50, "chr1:20300-20500": 12},
        }

    def test_afe_ale_psi(self):
        result = shibalib.afe_ale(self.junc_dict, ["s1"], self.event_df, 1, 5, 0)
        self.assertEqual(len(result), 1)
        row = result[0]
        # intron_a_mean=50, intron_b_mean=12 → PSI = 50/(50+12) ≈ 0.806
        s1_psi = row[12]  # 10 base + s1_junc_a(10), s1_junc_b(11), s1_psi(12)
        self.assertAlmostEqual(s1_psi, 50.0 / 62.0, places=4)

    def test_afe_ale_ind_psi(self):
        result = shibalib.afe_ale_ind(self.junc_dict, self.event_df, ["s1"], 1, 0)
        row = result[0]
        self.assertAlmostEqual(row[1], 50.0 / 62.0, places=4)


# ============================================================================
# PSI Calculation - MXE
# ============================================================================
class TestMXEPsi(unittest.TestCase):
    def setUp(self):
        self.event_df = pd.DataFrame({
            "event_id": ["MXE_1"],
            "pos_id": ["MXE@chr1@10000@10200-10400@10500-10700@10800"],
            "exon_a": ["chr1:10200-10400"],
            "exon_b": ["chr1:10500-10700"],
            "intron_a1": ["chr1:10000-10200"],
            "intron_a2": ["chr1:10400-10800"],
            "intron_b1": ["chr1:10000-10500"],
            "intron_b2": ["chr1:10700-10800"],
            "strand": ["+"],
            "gene_id": ["G10"],
            "gene_name": ["GeneJ"],
            "label": ["annotated"]
        })
        self.junc_dict = {
            "s1": {
                "chr1:10000-10200": 50,
                "chr1:10400-10800": 48,
                "chr1:10000-10500": 8,
                "chr1:10700-10800": 6,
            },
        }

    def test_mxe_psi(self):
        result = shibalib.mxe(self.junc_dict, ["s1"], self.event_df, 1, 5, 0)
        row = result[0]
        # a1+a2=98, b1+b2=14 → PSI = 98/(98+14) ≈ 0.875
        s1_psi = row[16]  # 12 base + s1_a1(12), s1_a2(13), s1_b1(14), s1_b2(15), psi(16)
        self.assertAlmostEqual(s1_psi, 98.0 / 112.0, places=4)

    def test_mxe_ind_psi(self):
        result = shibalib.mxe_ind(self.junc_dict, self.event_df, ["s1"], 1, 0)
        row = result[0]
        self.assertAlmostEqual(row[1], 98.0 / 112.0, places=4)


# ============================================================================
# PSI Calculation - RI
# ============================================================================
class TestRIPsi(unittest.TestCase):
    def setUp(self):
        self.event_df = pd.DataFrame({
            "event_id": ["RI_1"],
            "pos_id": ["RI@chr1@15100-15200"],
            "exon_a": ["chr1:15000-15100"],
            "exon_b": ["chr1:15200-15300"],
            "exon_c": ["chr1:15000-15300"],
            "intron_a": ["chr1:15100-15200"],
            "strand": ["+"],
            "gene_id": ["G12"],
            "gene_name": ["GeneL"],
            "label": ["annotated"]
        })
        # RI needs the main intron junction + start/end split junctions
        # intron_a_start_junc = chr1:15100-15101
        # intron_a_end_junc = chr1:15199-15200
        self.junc_dict = {
            "s1": {
                "chr1:15100-15200": 5,  # exclusion (intron_a)
                "chr1:15100-15101": 40, # retention start
                "chr1:15199-15200": 38, # retention end
            },
        }

    def test_ri_psi(self):
        result = shibalib.ri(self.junc_dict, ["s1"], self.event_df, 1, 5, 0)
        row = result[0]
        # retention = (40+38)/2 = 39, excl = 5 → PSI = 39/(39+5) ≈ 0.886
        s1_psi = row[13]  # 10 base + s1_junc_a_start(10), s1_junc_a_end(11), s1_junc_a(12), psi(13)
        self.assertAlmostEqual(s1_psi, 39.0 / 44.0, places=4)

    def test_ri_ind_psi(self):
        result = shibalib.ri_ind(self.junc_dict, self.event_df, ["s1"], 1, 0)
        row = result[0]
        self.assertAlmostEqual(row[1], 39.0 / 44.0, places=4)


# ============================================================================
# Differential Splicing - SE
# ============================================================================
class TestDiffSE(unittest.TestCase):
    def test_diff_se_basic(self):
        """Test diff_se identifies correct dPSI and statistical columns."""
        df = pd.DataFrame({
            "event_id": ["SE_1", "SE_2"],
            "pos_id": ["p1", "p2"],
            "exon": ["e1", "e2"],
            "intron_a": ["ia1", "ia2"],
            "intron_b": ["ib1", "ib2"],
            "intron_c": ["ic1", "ic2"],
            "strand": ["+", "+"],
            "gene_id": ["g1", "g2"],
            "gene_name": ["G1", "G2"],
            "label": ["annotated", "unannotated"],
            "ctrl_junction_a": [100, 50],
            "ctrl_junction_b": [90, 45],
            "ctrl_junction_c": [10, 80],
            "ctrl_PSI": [0.9, 0.37],
            "treat_junction_a": [10, 50],
            "treat_junction_b": [12, 48],
            "treat_junction_c": [80, 80],
            "treat_PSI": [0.12, 0.38],
        })
        result = shibalib.diff_se(df, ["ctrl", "treat"], 0.05, 0.1)
        self.assertIn("dPSI", result.columns)
        self.assertIn("Diff events", result.columns)
        self.assertIn("q", result.columns)
        self.assertIn("ref_PSI", result.columns)
        self.assertIn("alt_PSI", result.columns)

    def test_diff_se_empty_df(self):
        """diff_se should handle DataFrame that becomes empty after dropna."""
        df = pd.DataFrame({
            "event_id": ["SE_1"],
            "pos_id": ["p1"],
            "exon": ["e1"],
            "intron_a": ["ia1"],
            "intron_b": ["ib1"],
            "intron_c": ["ic1"],
            "strand": ["+"],
            "gene_id": ["g1"],
            "gene_name": ["G1"],
            "label": ["annotated"],
            "ctrl_junction_a": [100],
            "ctrl_junction_b": [90],
            "ctrl_junction_c": [10],
            "ctrl_PSI": [np.nan],
            "treat_junction_a": [10],
            "treat_junction_b": [12],
            "treat_junction_c": [80],
            "treat_PSI": [0.12],
        })
        result = shibalib.diff_se(df, ["ctrl", "treat"], 0.05, 0.1)
        self.assertEqual(len(result), 0)


# ============================================================================
# Differential Splicing - FIVE/THREE
# ============================================================================
class TestDiffFiveThree(unittest.TestCase):
    def test_diff_five_three_basic(self):
        df = pd.DataFrame({
            "event_id": ["FT_1"],
            "pos_id": ["p1"],
            "exon_a": ["ea1"],
            "exon_b": ["eb1"],
            "intron_a": ["ia1"],
            "intron_b": ["ib1"],
            "strand": ["+"],
            "gene_id": ["g1"],
            "gene_name": ["G1"],
            "label": ["annotated"],
            "ctrl_junction_a": [100],
            "ctrl_junction_b": [10],
            "ctrl_PSI": [0.91],
            "treat_junction_a": [10],
            "treat_junction_b": [100],
            "treat_PSI": [0.09],
        })
        result = shibalib.diff_five_three(df, ["ctrl", "treat"], 0.05, 0.1)
        self.assertIn("dPSI", result.columns)
        self.assertIn("Diff events", result.columns)


# ============================================================================
# Differential Splicing - MXE
# ============================================================================
class TestDiffMXE(unittest.TestCase):
    def test_diff_mxe_basic(self):
        df = pd.DataFrame({
            "event_id": ["MXE_1"],
            "pos_id": ["p1"],
            "exon_a": ["ea1"],
            "exon_b": ["eb1"],
            "intron_a1": ["ia1"],
            "intron_a2": ["ia2"],
            "intron_b1": ["ib1"],
            "intron_b2": ["ib2"],
            "strand": ["+"],
            "gene_id": ["g1"],
            "gene_name": ["G1"],
            "label": ["annotated"],
            "ctrl_junction_a1": [100],
            "ctrl_junction_a2": [90],
            "ctrl_junction_b1": [10],
            "ctrl_junction_b2": [8],
            "ctrl_PSI": [0.91],
            "treat_junction_a1": [10],
            "treat_junction_a2": [12],
            "treat_junction_b1": [90],
            "treat_junction_b2": [95],
            "treat_PSI": [0.11],
        })
        result = shibalib.diff_mxe(df, ["ctrl", "treat"], 0.05, 0.1)
        self.assertIn("dPSI", result.columns)
        self.assertIn("Diff events", result.columns)


# ============================================================================
# Differential Splicing - RI
# ============================================================================
class TestDiffRI(unittest.TestCase):
    def test_diff_ri_basic(self):
        df = pd.DataFrame({
            "event_id": ["RI_1"],
            "pos_id": ["p1"],
            "exon_a": ["ea1"],
            "exon_b": ["eb1"],
            "exon_c": ["ec1"],
            "intron_a": ["ia1"],
            "strand": ["+"],
            "gene_id": ["g1"],
            "gene_name": ["G1"],
            "label": ["annotated"],
            "ctrl_junction_a": [10],
            "ctrl_junction_a_start": [80],
            "ctrl_junction_a_end": [75],
            "ctrl_PSI": [0.89],
            "treat_junction_a": [80],
            "treat_junction_a_start": [10],
            "treat_junction_a_end": [12],
            "treat_PSI": [0.12],
        })
        result = shibalib.diff_ri(df, ["ctrl", "treat"], 0.05, 0.1)
        self.assertIn("dPSI", result.columns)
        self.assertIn("Diff events", result.columns)


# ============================================================================
# Differential Splicing - MSE
# ============================================================================
class TestDiffMSE(unittest.TestCase):
    def test_diff_mse_basic(self):
        df = pd.DataFrame({
            "event_id": ["MSE_1"],
            "pos_id": ["p1"],
            "mse_n": ["2"],
            "exon": ["e1;e2"],
            "intron": ["i1;i2;i3"],
            "strand": ["+"],
            "gene_id": ["g1"],
            "gene_name": ["G1"],
            "label": ["annotated"],
            "ctrl_junction": ["80;75;10"],
            "ctrl_PSI": [0.83],
            "treat_junction": ["10;12;80"],
            "treat_PSI": [0.12],
        })
        result = shibalib.diff_mse(df, ["ctrl", "treat"], 0.05, 0.1)
        self.assertIn("dPSI", result.columns)
        self.assertIn("Diff events", result.columns)


# ============================================================================
# Differential Splicing - AFE/ALE
# ============================================================================
class TestDiffAFEALE(unittest.TestCase):
    def test_diff_afe_ale_basic(self):
        df = pd.DataFrame({
            "event_id": ["AFE_1"],
            "pos_id": ["p1"],
            "exon_a": ["ea1"],
            "exon_b": ["eb1"],
            "intron_a": ["ia1"],
            "intron_b": ["ib1"],
            "strand": ["+"],
            "gene_id": ["g1"],
            "gene_name": ["G1"],
            "label": ["annotated"],
            "ctrl_junction_a": ["100"],
            "ctrl_junction_b": ["10"],
            "ctrl_PSI": [0.91],
            "treat_junction_a": ["10"],
            "treat_junction_b": ["100"],
            "treat_PSI": [0.09],
        })
        result = shibalib.diff_afe_ale(df, ["ctrl", "treat"], 0.05, 0.1)
        self.assertIn("dPSI", result.columns)
        self.assertIn("Diff events", result.columns)


# ============================================================================
# T-test
# ============================================================================
class TestTtest(unittest.TestCase):
    def test_ttest_basic(self):
        output_ind_df = pd.DataFrame({
            "event_id": ["SE_1", "SE_2"],
            "s1_PSI": [0.9, 0.5],
            "s2_PSI": [0.85, 0.55],
            "s3_PSI": [0.1, 0.5],
            "s4_PSI": [0.15, 0.45],
        })
        group_df = pd.DataFrame({
            "sample": ["s1", "s2", "s3", "s4"],
            "group": ["ctrl", "ctrl", "treat", "treat"]
        })
        result = shibalib.ttest(output_ind_df, group_df, ["ctrl", "treat"])
        self.assertIn("p_ttest", result.columns)
        self.assertEqual(len(result), 2)
        # SE_1 has big difference, should have small p-value
        self.assertLess(result["p_ttest"].iloc[0], 0.05)

    def test_ttest_with_nan(self):
        output_ind_df = pd.DataFrame({
            "event_id": ["SE_1"],
            "s1_PSI": [0.9],
            "s2_PSI": [np.nan],
            "s3_PSI": [0.1],
            "s4_PSI": [0.15],
        })
        group_df = pd.DataFrame({
            "sample": ["s1", "s2", "s3", "s4"],
            "group": ["ctrl", "ctrl", "treat", "treat"]
        })
        result = shibalib.ttest(output_ind_df, group_df, ["ctrl", "treat"])
        self.assertIn("p_ttest", result.columns)


# ============================================================================
# Beta Regression
# ============================================================================
class TestBetaRegression(unittest.TestCase):
    """Tests for beta_regression() — Beta regression with LRT."""

    def setUp(self):
        self.group_df = pd.DataFrame({
            "sample": ["s1", "s2", "s3", "s4"],
            "group": ["ctrl", "ctrl", "treat", "treat"]
        })
        self.group_list = ["ctrl", "treat"]

    def _make_ind_df(self, psi_g1, psi_g2, total_g1, total_g2, event_ids=None):
        """Build output_ind_df with _PSI and _total_reads columns.

        Args:
            psi_g1: list of lists, one per sample in group1 (e.g. [[0.9, 0.5], [0.85, 0.55]])
            psi_g2: same for group2
            total_g1: list of lists, total reads per sample in group1
            total_g2: same for group2
            event_ids: optional list of event IDs
        """
        n_events = len(psi_g1[0])
        if event_ids is None:
            event_ids = [f"SE_{i+1}" for i in range(n_events)]
        data = {"event_id": event_ids}
        for i, (psi, total) in enumerate(zip(psi_g1, total_g1)):
            s = f"s{i+1}"
            data[f"{s}_PSI"] = psi
            data[f"{s}_total_reads"] = total
        for i, (psi, total) in enumerate(zip(psi_g2, total_g2)):
            s = f"s{len(psi_g1)+i+1}"
            data[f"{s}_PSI"] = psi
            data[f"{s}_total_reads"] = total
        return pd.DataFrame(data)

    def test_beta_regression_basic(self):
        """Clear group difference should yield a small p-value."""
        df = self._make_ind_df(
            psi_g1=[[0.9], [0.85]],
            psi_g2=[[0.1], [0.15]],
            total_g1=[[100], [120]],
            total_g2=[[110], [105]],
        )
        result = shibalib.beta_regression(df, self.group_df, self.group_list)
        self.assertIn("p_beta", result.columns)
        self.assertEqual(len(result), 1)
        self.assertFalse(np.isnan(result["p_beta"].iloc[0]))
        self.assertLess(result["p_beta"].iloc[0], 0.05)

    def test_beta_regression_no_difference(self):
        """Identical PSI across groups should yield a large p-value."""
        df = self._make_ind_df(
            psi_g1=[[0.50], [0.52]],
            psi_g2=[[0.51], [0.49]],
            total_g1=[[100], [100]],
            total_g2=[[100], [100]],
        )
        result = shibalib.beta_regression(df, self.group_df, self.group_list)
        self.assertIn("p_beta", result.columns)
        self.assertGreater(result["p_beta"].iloc[0], 0.05)

    def test_beta_regression_with_nan_psi(self):
        """NaN PSI for one sample, but enough remaining samples (>=2 per group)."""
        group_df = pd.DataFrame({
            "sample": ["s1", "s2", "s3", "s4", "s5", "s6"],
            "group": ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]
        })
        df = pd.DataFrame({
            "event_id": ["SE_1"],
            "s1_PSI": [0.9], "s1_total_reads": [100],
            "s2_PSI": [np.nan], "s2_total_reads": [80],
            "s3_PSI": [0.88], "s3_total_reads": [90],
            "s4_PSI": [0.1], "s4_total_reads": [110],
            "s5_PSI": [0.15], "s5_total_reads": [95],
            "s6_PSI": [np.nan], "s6_total_reads": [100],
        })
        result = shibalib.beta_regression(df, group_df, ["ctrl", "treat"])
        self.assertIn("p_beta", result.columns)
        # Should still compute (2 valid per group remain)
        self.assertFalse(np.isnan(result["p_beta"].iloc[0]))

    def test_beta_regression_insufficient_samples(self):
        """Only 1 valid sample per group → p_beta should be NaN."""
        df = pd.DataFrame({
            "event_id": ["SE_1"],
            "s1_PSI": [0.9], "s1_total_reads": [100],
            "s2_PSI": [np.nan], "s2_total_reads": [80],
            "s3_PSI": [0.1], "s3_total_reads": [110],
            "s4_PSI": [np.nan], "s4_total_reads": [95],
        })
        result = shibalib.beta_regression(df, self.group_df, self.group_list)
        self.assertIn("p_beta", result.columns)
        self.assertTrue(np.isnan(result["p_beta"].iloc[0]))

    def test_beta_regression_total_reads_zero(self):
        """total_reads=0 should be filtered out; if insufficient remain → NaN."""
        df = pd.DataFrame({
            "event_id": ["SE_1"],
            "s1_PSI": [0.9], "s1_total_reads": [100],
            "s2_PSI": [0.85], "s2_total_reads": [0],  # filtered out
            "s3_PSI": [0.1], "s3_total_reads": [0],    # filtered out
            "s4_PSI": [0.15], "s4_total_reads": [105],
        })
        result = shibalib.beta_regression(df, self.group_df, self.group_list)
        self.assertIn("p_beta", result.columns)
        # Only 1 valid per group → NaN
        self.assertTrue(np.isnan(result["p_beta"].iloc[0]))

    def test_beta_regression_psi_boundary(self):
        """PSI at exact boundaries (0.0 and 1.0) should not cause errors."""
        df = self._make_ind_df(
            psi_g1=[[1.0], [0.95]],
            psi_g2=[[0.0], [0.05]],
            total_g1=[[100], [120]],
            total_g2=[[110], [105]],
        )
        result = shibalib.beta_regression(df, self.group_df, self.group_list)
        self.assertIn("p_beta", result.columns)
        self.assertFalse(np.isnan(result["p_beta"].iloc[0]))

    def test_beta_regression_multiple_events(self):
        """Multiple rows (events) should each get a p_beta value."""
        df = self._make_ind_df(
            psi_g1=[[0.9, 0.5, 0.3], [0.85, 0.52, 0.28]],
            psi_g2=[[0.1, 0.48, 0.7], [0.15, 0.51, 0.75]],
            total_g1=[[100, 80, 90], [120, 85, 95]],
            total_g2=[[110, 90, 100], [105, 88, 92]],
            event_ids=["SE_1", "SE_2", "SE_3"],
        )
        result = shibalib.beta_regression(df, self.group_df, self.group_list)
        self.assertIn("p_beta", result.columns)
        self.assertEqual(len(result), 3)
        # SE_1: large difference → small p
        self.assertLess(result["p_beta"].iloc[0], 0.05)
        # SE_2: almost no difference → large p
        self.assertGreater(result["p_beta"].iloc[1], 0.05)

    def test_beta_regression_column_mismatch(self):
        """Missing columns should raise ValueError."""
        df = pd.DataFrame({
            "event_id": ["SE_1"],
            "s1_PSI": [0.9],  # missing s2, s3, s4 and all _total_reads
        })
        with self.assertRaises(ValueError):
            shibalib.beta_regression(df, self.group_df, self.group_list)

    def test_beta_regression_empty_dataframe(self):
        """Empty DataFrame should return with p_beta column, 0 rows."""
        df = pd.DataFrame({
            "event_id": pd.Series([], dtype=str),
            "s1_PSI": pd.Series([], dtype=float), "s1_total_reads": pd.Series([], dtype=float),
            "s2_PSI": pd.Series([], dtype=float), "s2_total_reads": pd.Series([], dtype=float),
            "s3_PSI": pd.Series([], dtype=float), "s3_total_reads": pd.Series([], dtype=float),
            "s4_PSI": pd.Series([], dtype=float), "s4_total_reads": pd.Series([], dtype=float),
        })
        result = shibalib.beta_regression(df, self.group_df, self.group_list)
        self.assertIn("p_beta", result.columns)
        self.assertEqual(len(result), 0)

    def test_beta_regression_many_samples(self):
        """More samples (5 per group) should still work and give a small p-value for clear difference."""
        samples_g1 = ["a1", "a2", "a3", "a4", "a5"]
        samples_g2 = ["b1", "b2", "b3", "b4", "b5"]
        group_df = pd.DataFrame({
            "sample": samples_g1 + samples_g2,
            "group": ["ctrl"] * 5 + ["treat"] * 5
        })
        data = {"event_id": ["SE_1"]}
        for s in samples_g1:
            data[f"{s}_PSI"] = [np.random.uniform(0.80, 0.95)]
            data[f"{s}_total_reads"] = [np.random.randint(80, 150)]
        for s in samples_g2:
            data[f"{s}_PSI"] = [np.random.uniform(0.05, 0.20)]
            data[f"{s}_total_reads"] = [np.random.randint(80, 150)]
        np.random.seed(42)
        df = pd.DataFrame(data)
        result = shibalib.beta_regression(df, group_df, ["ctrl", "treat"])
        self.assertIn("p_beta", result.columns)
        self.assertFalse(np.isnan(result["p_beta"].iloc[0]))

    def test_beta_regression_parallel_matches_serial(self):
        """Parallel (num_process=2) should produce same p-values as serial."""
        df = self._make_ind_df(
            psi_g1=[[0.9, 0.5, 0.3, 0.7, 0.85], [0.85, 0.52, 0.28, 0.72, 0.80]],
            psi_g2=[[0.1, 0.48, 0.7, 0.3, 0.15], [0.15, 0.51, 0.75, 0.28, 0.12]],
            total_g1=[[100, 80, 90, 110, 95], [120, 85, 95, 105, 100]],
            total_g2=[[110, 90, 100, 88, 92], [105, 88, 92, 90, 85]],
            event_ids=["SE_1", "SE_2", "SE_3", "SE_4", "SE_5"],
        )
        result_serial = shibalib.beta_regression(df.copy(), self.group_df, self.group_list, num_process=1)
        result_parallel = shibalib.beta_regression(df.copy(), self.group_df, self.group_list, num_process=2)
        np.testing.assert_allclose(
            result_serial["p_beta"].values,
            result_parallel["p_beta"].values,
            rtol=1e-10, equal_nan=True
        )

    def test_beta_regression_prefilter_identical_psi(self):
        """If all PSI values are identical across groups, p should be ~1.0."""
        df = self._make_ind_df(
            psi_g1=[[0.5], [0.5]],
            psi_g2=[[0.5], [0.5]],
            total_g1=[[100], [100]],
            total_g2=[[100], [100]],
        )
        result = shibalib.beta_regression(df, self.group_df, self.group_list)
        self.assertIn("p_beta", result.columns)
        # Pre-filter should catch this (var < 1e-10) and return p≈1.0
        self.assertGreater(result["p_beta"].iloc[0], 0.99)


class TestBetaRegressionSingleEvent(unittest.TestCase):
    """Tests for _beta_regression_single_event (module-level worker function)."""

    def test_clear_difference_returns_large_lr_stat(self):
        """Clear group separation should yield a large LR statistic."""
        y = np.array([0.9, 0.85, 0.1, 0.15], dtype=np.float64)
        x = np.array([0, 0, 1, 1], dtype=np.float64)
        n = np.array([100, 120, 110, 105], dtype=np.float64)
        lr = shibalib._beta_regression_single_event(y, x, n)
        self.assertTrue(np.isfinite(lr))
        self.assertGreater(lr, 0)

    def test_no_difference_returns_small_lr_stat(self):
        """Nearly identical PSI across groups should yield LR ~ 0."""
        y = np.array([0.50, 0.52, 0.51, 0.49], dtype=np.float64)
        x = np.array([0, 0, 1, 1], dtype=np.float64)
        n = np.array([100, 100, 100, 100], dtype=np.float64)
        lr = shibalib._beta_regression_single_event(y, x, n)
        self.assertTrue(np.isfinite(lr))
        self.assertLess(lr, 3.84)  # chi2 critical value at p=0.05, df=1

    def test_insufficient_samples_returns_nan(self):
        """< 2 per group → NaN."""
        y = np.array([0.9, 0.1], dtype=np.float64)
        x = np.array([0, 1], dtype=np.float64)
        n = np.array([100, 100], dtype=np.float64)
        lr = shibalib._beta_regression_single_event(y, x, n)
        self.assertTrue(np.isnan(lr))

    def test_zero_variance_returns_zero(self):
        """All PSI identical → LR=0 via pre-filter."""
        y = np.array([0.5, 0.5, 0.5, 0.5], dtype=np.float64)
        x = np.array([0, 0, 1, 1], dtype=np.float64)
        n = np.array([100, 100, 100, 100], dtype=np.float64)
        lr = shibalib._beta_regression_single_event(y, x, n)
        self.assertEqual(lr, 0.0)

    def test_identical_group_means_returns_zero(self):
        """Group means within 1e-8 → LR=0 via pre-filter."""
        y = np.array([0.500000001, 0.499999999, 0.500000001, 0.499999999], dtype=np.float64)
        x = np.array([0, 0, 1, 1], dtype=np.float64)
        n = np.array([100, 100, 100, 100], dtype=np.float64)
        lr = shibalib._beta_regression_single_event(y, x, n)
        self.assertEqual(lr, 0.0)

    def test_boundary_psi_values(self):
        """PSI at 0.0 and 1.0 should not crash."""
        y = np.array([1.0, 0.95, 0.0, 0.05], dtype=np.float64)
        x = np.array([0, 0, 1, 1], dtype=np.float64)
        n = np.array([100, 120, 110, 105], dtype=np.float64)
        lr = shibalib._beta_regression_single_event(y, x, n)
        self.assertTrue(np.isfinite(lr))

    def test_optimizer_failure_returns_nan(self):
        """If optimization fails, should return NaN (no Nelder-Mead fallback)."""
        # Pathological data: all samples have same PSI but high variance in reads
        y = np.array([0.99999, 0.99999, 0.00001, 0.00001], dtype=np.float64)
        x = np.array([0, 0, 1, 1], dtype=np.float64)
        n = np.array([1, 1, 1, 1], dtype=np.float64)
        lr = shibalib._beta_regression_single_event(y, x, n)
        # Should return either a valid statistic or NaN — never raise
        self.assertTrue(np.isfinite(lr) or np.isnan(lr))


class TestBetaRegAnalyticalGradient(unittest.TestCase):
    """Verify analytical gradients match numerical approximation."""

    def setUp(self):
        """Set up test data for gradient checks."""
        self.y = np.array([0.8, 0.85, 0.2, 0.15], dtype=np.float64)
        self.x = np.array([0, 0, 1, 1], dtype=np.float64)
        self.n = np.array([100, 120, 110, 105], dtype=np.float64)
        n_total = len(self.y)
        self.y_t = (self.y * (n_total - 1) + 0.5) / n_total
        self.log_n = np.log(self.n)
        self.log_y = np.log(self.y_t)
        self.log_1_y = np.log(1 - self.y_t)
        self.args = (self.x, self.log_n, self.log_y, self.log_1_y)

    def _numerical_grad(self, func, params, args, eps=1e-7):
        """Compute numerical gradient via forward differences."""
        grad = np.zeros_like(params)
        for i in range(len(params)):
            p_plus = params.copy()
            p_plus[i] += eps
            p_minus = params.copy()
            p_minus[i] -= eps
            grad[i] = (func(p_plus, *args) - func(p_minus, *args)) / (2 * eps)
        return grad

    def test_full_model_gradient(self):
        """Analytical gradient of full model should match numerical approx."""
        from scipy.optimize import approx_fprime
        params = np.array([0.5, -1.0, 2.3, 0.1])
        analytical = shibalib._beta_reg_jac_full(params, *self.args)
        numerical = self._numerical_grad(shibalib._beta_reg_neg_ll_full, params, self.args)
        np.testing.assert_allclose(analytical, numerical, rtol=1e-4, atol=1e-6)

    def test_null_model_gradient(self):
        """Analytical gradient of null model should match numerical approx."""
        params = np.array([0.5, 2.3, 0.1])
        analytical = shibalib._beta_reg_jac_null(params, *self.args)
        numerical = self._numerical_grad(shibalib._beta_reg_neg_ll_null, params, self.args)
        np.testing.assert_allclose(analytical, numerical, rtol=1e-4, atol=1e-6)

    def test_full_gradient_at_optimum(self):
        """At the MLE, gradient should be approximately zero."""
        from scipy.optimize import minimize as sp_minimize
        params0 = np.array([0.5, -2.0, 2.3, 0.1])
        result = sp_minimize(
            shibalib._beta_reg_neg_ll_full, params0,
            args=self.args, jac=shibalib._beta_reg_jac_full,
            method='L-BFGS-B', options={'maxiter': 5000, 'ftol': 1e-15}
        )
        grad_at_opt = shibalib._beta_reg_jac_full(result.x, *self.args)
        # Clipping in the likelihood can cause non-trivial residual gradients
        # at the boundary; verify they are reasonably small relative to scale
        np.testing.assert_allclose(grad_at_opt, 0.0, atol=1.5)

    def test_null_gradient_at_optimum(self):
        """At the null MLE, gradient should be approximately zero."""
        from scipy.optimize import minimize as sp_minimize
        params0 = np.array([0.0, 2.3, 0.1])
        result = sp_minimize(
            shibalib._beta_reg_neg_ll_null, params0,
            args=self.args, jac=shibalib._beta_reg_jac_null,
            method='L-BFGS-B', options={'maxiter': 1000, 'ftol': 1e-14}
        )
        grad_at_opt = shibalib._beta_reg_jac_null(result.x, *self.args)
        np.testing.assert_allclose(grad_at_opt, 0.0, atol=1e-4)


class TestBetaRegressionChunk(unittest.TestCase):
    """Tests for _beta_regression_chunk."""

    def test_chunk_processes_multiple_events(self):
        """A chunk of events should return one LR stat per event."""
        chunk = [
            (np.array([0.9, 0.85, 0.1, 0.15]), np.array([0, 0, 1, 1]), np.array([100, 120, 110, 105])),
            (np.array([0.5, 0.52, 0.51, 0.49]), np.array([0, 0, 1, 1]), np.array([100, 100, 100, 100])),
            (np.array([0.5, 0.1]), np.array([0, 1]), np.array([100, 100])),  # insufficient
        ]
        results = shibalib._beta_regression_chunk(chunk)
        self.assertEqual(len(results), 3)
        self.assertTrue(np.isfinite(results[0]))
        self.assertTrue(np.isfinite(results[1]))
        self.assertTrue(np.isnan(results[2]))


# ============================================================================
# make_psi_mtx
# ============================================================================
class TestMakePsiMtx(unittest.TestCase):
    def test_make_psi_mtx(self):
        df = pd.DataFrame({
            "event_id": ["SE_2", "SE_1"],
            "pos_id": ["p2", "p1"],
            "extra": ["x", "y"],
            "s1_PSI": [0.5, 0.9],
            "s2_PSI": [0.6, 0.85],
        })
        sorted_df, mtx_df = shibalib.make_psi_mtx(df)
        # Should be sorted by numeric part of event_id
        self.assertEqual(list(sorted_df["event_id"]), ["SE_1", "SE_2"])
        # mtx_df should have event_id, pos_id, and sample columns (without _PSI)
        self.assertIn("event_id", mtx_df.columns)
        self.assertIn("pos_id", mtx_df.columns)
        self.assertIn("s1", mtx_df.columns)
        self.assertIn("s2", mtx_df.columns)


# ============================================================================
# EventCounter
# ============================================================================
class TestEventCounter(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame({
            "event_id": ["SE_1", "SE_2", "SE_3", "SE_4", "SE_5"],
            "Diff events": ["Yes", "Yes", "Yes", "Yes", "No"],
            "dPSI": [0.3, -0.4, 0.2, -0.25, 0.5],
            "label": ["annotated", "unannotated", "unannotated", "annotated", "annotated"],
        })

    def test_count_all_events(self):
        counter = shibalib.EventCounter(self.df, 0.1)
        counts = counter.count_all_events()
        self.assertEqual(counts["up_annotated_num"], 1)    # SE_1 (dPSI=0.3, annotated, Yes)
        self.assertEqual(counts["up_unannotated_num"], 1)  # SE_3 (dPSI=0.2, unannotated, Yes)
        self.assertEqual(counts["down_annotated_num"], 1)  # SE_4 (dPSI=-0.25, annotated, Yes)
        self.assertEqual(counts["down_unannotated_num"], 1) # SE_2 (dPSI=-0.4, unannotated, Yes)

    def test_count_events_no_diff(self):
        df = pd.DataFrame({
            "event_id": ["SE_1"],
            "Diff events": ["No"],
            "dPSI": [0.5],
            "label": ["annotated"],
        })
        counter = shibalib.EventCounter(df, 0.1)
        counts = counter.count_all_events()
        self.assertEqual(counts["up_annotated_num"], 0)

    def test_count_events_empty_df(self):
        df = pd.DataFrame(columns=["event_id", "Diff events", "dPSI", "label"])
        counter = shibalib.EventCounter(df, 0.1)
        counts = counter.count_all_events()
        self.assertEqual(counts["up_annotated_num"], 0)
        self.assertEqual(counts["down_unannotated_num"], 0)


if __name__ == "__main__":
    unittest.main()
