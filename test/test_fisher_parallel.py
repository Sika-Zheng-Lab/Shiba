"""Parallel Fisher family correction and observed-PSI regression checks."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from multiprocessing.shared_memory import SharedMemory

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lib import factorial, shibalib

DATA = Path(__file__).parent / "data"
FUNCTIONS = {
    "SE": (shibalib.diff_se, shibalib.col_se),
    "MSE": (shibalib.diff_mse, shibalib.col_mse),
    "MXE": (shibalib.diff_mxe, shibalib.col_mxe),
    "RI": (shibalib.diff_ri, shibalib.col_ri),
    "FIVE": (shibalib.diff_five_three, shibalib.col_five_three_afe_ale),
    "THREE": (shibalib.diff_five_three, shibalib.col_five_three_afe_ale),
    "AFE": (shibalib.diff_afe_ale, shibalib.col_five_three_afe_ale),
    "ALE": (shibalib.diff_afe_ale, shibalib.col_five_three_afe_ale),
}


def fisher_fixture(kind, n=384):
    events = pd.read_csv(DATA / f"EVENT_{kind}.txt", sep="\t")
    df = events.iloc[np.arange(n) % len(events)].reset_index(drop=True).copy()
    df["event_id"] = [f"{kind}_{i}" for i in range(n)]
    rng = np.random.default_rng(2026)
    for col in FUNCTIONS[kind][1](["ctrl", "treat"], False):
        if col in df:
            continue
        if col.endswith("_PSI"):
            df[col] = rng.uniform(0, 1, n)
        else:
            counts = rng.integers(0, 60, (n, 3))
            if kind == "MSE":
                df[col] = [";".join(map(str, row)) for row in counts]
            elif kind in ("AFE", "ALE"):
                df[col] = [";".join(map(str, row[:1 + i % 3])) for i, row in enumerate(counts)]
            else:
                df[col] = counts[:, 0]
    return df


class TestParallelFisher(unittest.TestCase):
    def test_all_event_types_and_global_fdr(self):
        for kind, (func, _) in FUNCTIONS.items():
            with self.subTest(kind=kind):
                df = fisher_fixture(kind)
                # An entirely invalid task plus isolated missing rows.
                df.loc[:63, "ctrl_PSI"] = np.nan
                df.loc[200, "treat_PSI"] = np.nan
                serial = func(df, ["ctrl", "treat"], .05, .1)
                parallel = shibalib.fisher_event_tests(df, ["ctrl", "treat"], func, 2, .05, .1)
                pd.testing.assert_frame_equal(serial, parallel, check_exact=True)
                p = parallel["p" if kind in ("FIVE", "THREE") else "p_maximum"]
                q = multipletests(p, method="fdr_bh")[1]
                q[q == 0] = 1e-323
                np.testing.assert_array_equal(parallel.q, q)
                self.assertTrue((parallel["Diff events"] == "Yes").any())
                self.assertTrue((parallel["Diff events"] == "No").any())

    def test_empty_and_all_missing(self):
        for kind, (func, _) in FUNCTIONS.items():
            df = fisher_fixture(kind)
            df["ctrl_PSI"] = np.nan
            for frame in (df, df.iloc[:0]):
                with self.subTest(kind=kind, rows=len(frame)):
                    pd.testing.assert_frame_equal(
                        func(frame, ["ctrl", "treat"], .05, .1),
                        shibalib.fisher_event_tests(frame, ["ctrl", "treat"], func, 2, .05, .1))

    def test_cache_matches_scipy_including_boundaries(self):
        shibalib._cached_fisher_exact.cache_clear()
        for table in ([[0, 0], [0, 0]], [[20, 0], [0, 20]],
                      [[0, 10], [5, 20]], [[11, 8], [9, 3]]):
            expected = tuple(fisher_exact(table, alternative="two-sided"))
            np.testing.assert_equal(shibalib._fisher_exact(table), expected)
            np.testing.assert_equal(shibalib._fisher_exact(table), expected)
        self.assertEqual(shibalib._cached_fisher_exact.cache_info().hits, 4)

    def test_diff_event_uses_parallel_dispatch(self):
        df = fisher_fixture("SE")
        expected = shibalib.diff_se(df, ["ctrl", "treat"], .05, .1).sort_values(["Diff events", "q"], ascending=[False, True]).reset_index(drop=True)
        result = shibalib.diff_event(None, df, None, None, ["ctrl", "treat"], [],
                                    shibalib.diff_se, None, 2, .05, .1, False, False)
        pd.testing.assert_frame_equal(result.reset_index(drop=True), expected)


class TestObservedPSIParallel(unittest.TestCase):
    def setUp(self):
        self.events = shibalib.read_events(str(DATA))
        self.junctions = shibalib.JunctionData.from_dataframe(pd.read_csv(DATA / "junctions.bed", sep="\t"))
        self.samples = self.junctions.sample_ids

    def test_serial_parallel_identical_all_types(self):
        serial_matrix, serial_tables = factorial.observed_psi(self.events, self.junctions, self.samples, 10, 1)
        matrix, tables = factorial.observed_psi(self.events, self.junctions, self.samples, 10, 2)
        pd.testing.assert_frame_equal(matrix, serial_matrix, check_exact=True)
        for kind in self.events:
            pd.testing.assert_frame_equal(tables[kind], serial_tables[kind], check_exact=True)

    def test_empty_events(self):
        events = {kind: frame.iloc[:0] for kind, frame in self.events.items()}
        serial, _ = factorial.observed_psi(events, self.junctions, self.samples, 10, 1)
        parallel, _ = factorial.observed_psi(events, self.junctions, self.samples, 10, 2)
        pd.testing.assert_frame_equal(serial, parallel)

    def test_shared_memory_released_on_failure(self):
        shm, info = self.junctions.to_shared_memory()
        with patch.object(self.junctions, "to_shared_memory", return_value=(shm, info)), \
             patch.object(shibalib, "make_psi_table_sample", side_effect=RuntimeError("worker failed")):
            with self.assertRaisesRegex(RuntimeError, "worker failed"):
                factorial.observed_psi(self.events, self.junctions, self.samples, 10, 2)
        with self.assertRaises(FileNotFoundError):
            SharedMemory(name=info["shm_name"])


if __name__ == "__main__":
    unittest.main()
