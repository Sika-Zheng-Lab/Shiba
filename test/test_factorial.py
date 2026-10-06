"""Regression contract tests, with optional real glmmTMB integration tests."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lib import factorial, shibalib
from lib.event_components import component_pairs, extract_components, EVENT_TYPES

ROOT = Path(__file__).resolve().parents[1]


def has_backend():
    if not shutil.which("Rscript"):
        return False
    result = subprocess.run(["Rscript", "--vanilla", "-e",
                             'quit(status=if(requireNamespace("glmmTMB",quietly=TRUE) && requireNamespace("jsonlite",quietly=TRUE)) 0 else 1)'],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return result.returncode == 0


def make_fixture(directory, replicates=12):
    """Shared raw junctions; BB samples, including null and contradictory SEs."""
    directory = Path(directory)
    (directory / "events").mkdir()
    for kind in EVENT_TYPES:
        frame = pd.read_csv(ROOT / "test" / "data" / f"EVENT_{kind}.txt", sep="\t")
        frame.iloc[:0].to_csv(directory / "events" / f"EVENT_{kind}.txt", sep="\t", index=False)
    meta = pd.DataFrame([(f"s{g}{t}{r}", g, t, r % 2, r / replicates)
                         for g in ["WT", "KO"] for t in ["Control", "Drug"] for r in range(replicates)],
                        columns=["sample", "genotype", "treatment", "batch", "time"])
    meta.to_csv(directory / "samples.tsv", sep="\t", index=False)
    rng = np.random.default_rng(716)
    events, junctions = [], []
    pattern = ((meta.genotype == "KO") & (meta.treatment == "Drug")).to_numpy()
    for i, name in enumerate(["interaction", "null", "opposite", "missing", "boundary"], 1):
        base = 1000 * i
        a, b, c = f"chr1:{base}-{base+100}", f"chr1:{base+200}-{base+300}", f"chr1:{base}-{base+300}"
        events.append([f"SE_{i}", f"SE@{i}", f"chr1:{base+100}-{base+200}", a, b, c, "+", name, name, "annotated"])
        mu = np.where(pattern, .8, .2) if name != "null" else np.full(len(meta), .35)
        p = rng.beta(mu * 40, (1-mu) * 40)
        inc = rng.binomial(200, p); exc = 200 - inc
        inc2 = inc.copy()
        if name == "opposite":
            opposite = rng.beta(np.where(pattern, .05, .8) * 40, (1 - np.where(pattern, .05, .8)) * 40)
            inc2 = rng.poisson(exc * opposite / (1 - opposite))
        if name == "missing":
            inc[:] = 0; inc2[:] = 0; exc[:] = 0
        if name == "boundary":
            inc[:] = 0; inc2[:] = 0
        for key, values in [(a, inc), (b, inc2), (c, exc)]:
            chrom, coords = key.split(":"); start, end = coords.split("-")
            junctions.append([chrom, start, end, key, *values])
    columns = ["event_id", "pos_id", "exon", "intron_a", "intron_b", "intron_c", "strand", "gene_id", "gene_name", "label"]
    pd.DataFrame(events, columns=columns).to_csv(directory / "events" / "EVENT_SE.txt", sep="\t", index=False)
    pd.DataFrame(junctions, columns=["chr", "start", "end", "ID", *meta['sample']]).to_csv(directory / "junctions.bed", sep="\t", index=False)
    return meta


class TestComponents(unittest.TestCase):
    def test_all_event_definitions(self):
        expectations = {"SE": 2, "FIVE": 1, "THREE": 1, "MXE": 4, "RI": 2, "MSE": 3, "AFE": 1, "ALE": 1}
        for kind, count in expectations.items():
            row = pd.read_csv(ROOT / "test" / "data" / f"EVENT_{kind}.txt", sep="\t").iloc[0]
            self.assertEqual(len(component_pairs(kind, row)), count, kind)
        row = {"event_id": "AFE_1", "intron_a": "a1;a2", "intron_b": "b1;b2;b3"}
        self.assertEqual(component_pairs("AFE", row), [(a, b) for a in ["a1", "a2"] for b in ["b1", "b2", "b3"]])

    def test_integer_counts_and_shared_exclusion(self):
        events = pd.DataFrame([dict(event_id="SE_1", intron_a="a", intron_b="b", intron_c="c")])
        jd = shibalib.JunctionData(np.array([[3, 8], [11, 4]]), {"a": 0, "c": 1}, ["s1", "s2"])
        counts, definitions = extract_components("SE", events, jd, ["s2", "s1"])
        self.assertEqual(counts.inclusion.tolist(), [8, 3, 0, 0])
        self.assertEqual(counts.exclusion.tolist(), [4, 11, 4, 11])
        self.assertEqual(definitions.exclusion_junction.tolist(), ["c", "c"])
        self.assertTrue(pd.api.types.is_integer_dtype(counts.inclusion))

    def test_ri_retention_orientation(self):
        row = dict(event_id="RI_1", intron_a="chr1:10-30")
        self.assertEqual(component_pairs("RI", row), [("chr1:10-11", "chr1:10-30"), ("chr1:29-30", "chr1:10-30")])

    def test_invalid_components(self):
        with self.assertRaises(ValueError):
            component_pairs("SE", dict(event_id="SE_1", intron_a="a", intron_b="a", intron_c="c"))


class TestCoverageFilter(unittest.TestCase):
    def fixture(self):
        samples = ["a1", "a2", "b1", "b2", "b3"]
        groups = {"a": samples[:2], "b": samples[2:]}
        # Exactly half in a, and ceil(3/2) in b; 9 reads is below threshold.
        totals = {"pass": [10, 9, 10, 10, 9], "fail": [10, 10, 10, 9, 9]}
        rows = [["SE", event, component, sample, total // 2, total - total // 2]
                for event, values in totals.items() for component in ["c1", "c2"]
                for sample, total in zip(samples, values)]
        return pd.DataFrame(rows, columns=[*factorial.KEYS, "component_id", "sample", "inclusion", "exclusion"]), groups

    def test_threshold_each_group_and_odd_group_size(self):
        counts, groups = self.fixture()
        result = factorial.filter_events(counts, groups)
        self.assertEqual(result.groupby("event_id").filter_pass.first().to_dict(), {"fail": False, "pass": True})
        self.assertEqual(result[result.filter_group == "b"].n_required.tolist(), [2, 2])
        self.assertEqual(result[result.event_id == "pass"].n_eligible.tolist(), [1, 2])
        self.assertFalse(factorial.filter_events(counts, groups, fraction=1).filter_pass.any())
        self.assertTrue(factorial.filter_events(counts, groups, minimum_reads=9).filter_pass.all())

    def test_same_samples_must_cover_all_components(self):
        counts, _ = self.fixture()
        counts = counts[(counts.event_id == "pass") & counts['sample'].isin(["a1", "a2"])].copy()
        counts["exclusion"] = 0
        counts["inclusion"] = [10, 0, 0, 10]
        result = factorial.filter_events(counts, {"a": ["a1", "a2"]})
        self.assertFalse(result.filter_pass.any())
        self.assertEqual(result.n_eligible.iloc[0], 0)

    def test_group_resolution(self):
        meta = pd.DataFrame(dict(sample=["s1", "s2", "s3", "s4"], group=["a", "a", "b", "b"],
                                 genotype=["WT", "WT", "KO", "KO"], treatment=["C", "D", "C", "D"]))
        design = {"levels": {"genotype": ["WT", "KO"], "treatment": ["C", "D"]}}
        columns, groups = factorial.resolve_filter_groups(meta, design, [])
        self.assertEqual(columns, ["genotype", "treatment"])
        self.assertEqual(len(groups), 4)  # An unrelated group column does not override design cells.
        columns, groups = factorial.resolve_filter_groups(meta.drop(columns="group"), design, [])
        self.assertEqual(columns, ["genotype", "treatment"])
        self.assertEqual(len(groups), 4)
        columns, groups = factorial.resolve_filter_groups(meta, design, ["treatment"])
        self.assertEqual(list(groups.values()), [["s1", "s3"], ["s2", "s4"]])
        columns, groups = factorial.resolve_filter_groups(meta, design, ["group"])
        self.assertEqual(list(groups.values()), [["s1", "s2"], ["s3", "s4"]])
        columns, groups = factorial.resolve_filter_groups(meta, {}, [])
        self.assertEqual(list(groups.values()), [meta['sample'].tolist()])
        for columns in [["missing"], ["sample"], ["group", "group"]]:
            with self.assertRaises(ValueError):
                factorial.resolve_filter_groups(meta, design, columns)
        meta.loc[0, "group"] = ""
        with self.assertRaises(ValueError):
            factorial.resolve_filter_groups(meta, design, ["group"])

    def test_filtered_events_excluded_from_bh_and_by(self):
        counts, groups = self.fixture()
        filtering = factorial.filter_events(counts, groups)
        definitions = counts[factorial.KEYS + ["component_id"]].drop_duplicates()
        stats = definitions[definitions.event_id == "pass"].copy()
        stats["contrast"] = "test"; stats["status"] = "ok"
        stats["estimate"] = 1.; stats["p_component"] = .04; stats["n_samples"] = 3
        for adjustment in ["BH", "BY"]:
            result = factorial.aggregate_events(stats, definitions, ["test"], adjust=adjustment,
                                                filtering=filtering).set_index("event_id")
            self.assertEqual(result.loc["pass", "q_event"], .04)
            self.assertEqual(result.loc["fail", "status"], "filtered_low_coverage")
            self.assertTrue(np.isnan(result.loc["fail", "q_event"]))
            self.assertTrue(np.isnan(result.loc["fail", "p_event"]))

    def test_fraction_validation_and_defaults(self):
        from lib.general import validate_config_types
        parser = argparse.ArgumentParser()
        factorial.add_arguments(parser)
        args = parser.parse_args([])
        self.assertEqual(args.min_sample_fraction, .5)
        args.formula = "group"; args.sample_metadata = "samples.tsv"
        args.onlypsi = args.onlypsi_group = args.beta_binomial = args.ttest = False
        args.coef = ["groupB"]; args.num_process = 1; args.minimum_reads = 10; args.fdr = .05
        for fraction in [0, -0.1, 1.01, float("nan"), float("inf")]:
            args.min_sample_fraction = fraction
            with self.assertRaisesRegex(ValueError, "min-sample-fraction"):
                factorial.validate_options(args)
            self.assertTrue(validate_config_types(dict(stat_method="beta-binomial", min_sample_fraction=fraction)))
        args.min_sample_fraction = .5
        factorial.validate_options(args)


class TestAggregation(unittest.TestCase):
    def test_pipeline_config_without_pairwise_groups(self):
        from lib.general import validate_config
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "experiment.tsv").write_text("sample\tbam\tgroup\na\t/path/a.bam\tall\n")
            cfg = dict(stat_method="beta-binomial", formula="genotype * treatment",
                       coef=["genotypeKO:treatmentDrug"], experiment_table=str(tmp / "experiment.tsv"))
            self.assertEqual(validate_config(cfg), [])
            cfg["filter_group"] = [123]
            self.assertTrue(any("filter_group" in e for e in validate_config(cfg)))
            cfg["filter_group"] = ["genotype", "treatment"]
            self.assertEqual(validate_config(cfg), [])
            cfg["formula"] = ""
            self.assertTrue(any("requires formula" in e for e in validate_config(cfg)))

    def test_duplicate_column_names_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.tsv"
            path.write_text("sample\tcondition\tcondition\na\tWT\tKO\n")
            with self.assertRaises(ValueError):
                factorial.read_metadata(path, ["a"])

    def fixture(self):
        definitions = pd.DataFrame([(kind, event, c) for kind, event in [("SE", "SE_1"), ("MXE", "MXE_1")]
                                    for c in ["c1", "c2"]], columns=factorial.KEYS + ["component_id"])
        stats = definitions.copy()
        stats["contrast"] = "interaction"; stats["status"] = "ok"
        stats["estimate"] = [1., 2., 3., 4.]; stats["p_component"] = [.001, .02, .003, .005]
        stats["n_samples"] = 8
        return stats, definitions

    def test_maximum_and_fdr_across_types(self):
        stats, definitions = self.fixture()
        result = factorial.aggregate_events(stats, definitions, ["interaction"])
        np.testing.assert_allclose(result.p_event, [.02, .005])
        np.testing.assert_allclose(result.q_event, [.02, .01])
        self.assertTrue(result["Diff events"].eq("Yes").all())

    def test_opposite_or_missing_never_pass(self):
        stats, definitions = self.fixture()
        stats.loc[1, "estimate"] = -2
        stats.loc[3, "status"] = "null_fit_failed"
        result = factorial.aggregate_events(stats, definitions, ["interaction"])
        self.assertEqual(result.status.tolist(), ["inconsistent_direction", "untestable"])
        self.assertEqual(result.p_event.iloc[0], 1)
        self.assertTrue(np.isnan(result.p_event.iloc[1]))
        self.assertTrue(result["Diff events"].eq("No").all())
        result = factorial.aggregate_events(stats.iloc[:3], definitions, ["interaction"])
        self.assertEqual(result.status.iloc[1], "untestable")

    def test_missing_predictions_not_averaged_away(self):
        _, definitions = self.fixture()
        predictions = pd.DataFrame([["SE", "SE_1", "c1", "profile", "g1", .8]], columns=factorial.PRED_COLUMNS)
        result = factorial.aggregate_predictions(predictions, definitions)
        self.assertEqual(result.status.iloc[0], "incomplete_components")
        self.assertTrue(np.isnan(result.value.iloc[0]))

    def test_failures_remain_in_fdr_family(self):
        stats, definitions = self.fixture()
        stats.loc[2:, "status"] = "full_fit_failed"
        result = factorial.aggregate_events(stats, definitions, ["interaction"])
        self.assertAlmostEqual(result.q_event.iloc[0], .04)

    def test_metadata_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meta.tsv"
            path.write_text("sample\tcondition\nb\tKO\na\tWT\n")
            self.assertEqual(factorial.read_metadata(path, ["a", "b"]).condition.tolist(), ["WT", "KO"])
            with self.assertRaises(ValueError):
                factorial.read_metadata(path, ["a"])

    def test_config_arguments(self):
        args = factorial.config_arguments(dict(stat_method="beta-binomial", formula="a * b", coef="aKO:bDrug",
                                               reference_levels={"a": "WT"}, filter_group=["a", "b"], min_sample_fraction=.75))
        self.assertIn("a * b", args)
        self.assertEqual(args[args.index("--coef") + 1], "aKO:bDrug")
        self.assertEqual(args.count("--filter-group"), 2)
        self.assertEqual(args[args.index("--min-sample-fraction") + 1], "0.75")
        self.assertEqual(factorial.config_arguments({}), [])


@unittest.skipUnless(has_backend(), "R packages glmmTMB/jsonlite unavailable")
class TestGlmmTMB(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.metadata = make_fixture(cls.root)
        cls.base = [sys.executable, str(ROOT / "src" / "psi.py"), str(cls.root / "junctions.bed"),
                    str(cls.root / "events"), str(cls.root / "out"), "--stat-method", "beta-binomial",
                    "--sample-metadata", str(cls.root / "samples.tsv"), "--formula", "genotype * treatment + batch",
                    "--reference-level", "genotype=WT", "--reference-level", "treatment=Control",
                    "--coef", "genotypeKO:treatmentDrug"]
        result = subprocess.run(cls.base, capture_output=True, text=True)
        if result.returncode:
            raise AssertionError(result.stderr)
        cls.events = pd.read_csv(cls.root / "out" / "event_statistics.tsv", sep="\t")

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_known_interaction_and_consistency(self):
        frame = self.events.set_index("event_id")
        self.assertEqual(frame.loc["SE_1", "Diff events"], "Yes")
        self.assertGreater(frame.loc["SE_1", "estimate_min"], 0)
        self.assertEqual(frame.loc["SE_2", "Diff events"], "No")
        self.assertEqual(frame.loc["SE_3", "status"], "inconsistent_direction")
        self.assertEqual(frame.loc["SE_4", "status"], "filtered_low_coverage")
        self.assertEqual(frame.loc["SE_5", "status"], "untestable")
        prediction = pd.read_csv(self.root / "out" / "model_predictions.tsv", sep="\t")
        dd = prediction[(prediction.event_id == "SE_1") & prediction.name.str.startswith("delta_delta:")]
        self.assertEqual(len(dd), 1)
        self.assertGreater(dd.value.iloc[0], .4)

    def test_matches_direct_r_full_and_reduced_formula(self):
        script = self.root / "reference.R"
        script.write_text('''
suppressPackageStartupMessages(library(glmmTMB))
a <- commandArgs(TRUE)
d <- read.delim(a[1], check.names=FALSE)
j <- read.delim(a[2], check.names=FALSE)
d$genotype <- relevel(factor(d$genotype), "WT")
d$treatment <- relevel(factor(d$treatment), "Control")
d$k <- as.numeric(j[1, d$sample]); d$e <- as.numeric(j[3, d$sample])
full <- glmmTMB(cbind(k,e) ~ genotype*treatment+batch, data=d, family=betabinomial())
null <- glmmTMB(cbind(k,e) ~ genotype+treatment+batch, data=d, family=betabinomial())
write.table(data.frame(estimate=fixef(full)$cond["genotypeKO:treatmentDrug"],
    full=as.numeric(logLik(full)), null=as.numeric(logLik(null))), a[3],
    row.names=FALSE, sep="\\t")
''')
        reference = self.root / "reference.tsv"
        run = subprocess.run(["Rscript", "--vanilla", str(script), str(self.root / "samples.tsv"),
                              str(self.root / "junctions.bed"), str(reference)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        ref = pd.read_csv(reference, sep="\t").iloc[0]
        components = pd.read_csv(self.root / "out" / "component_statistics.tsv", sep="\t")
        result = components[components.event_id == "SE_1"].iloc[0]
        self.assertAlmostEqual(result.estimate, ref.estimate, places=4)
        self.assertAlmostEqual(result.loglik_full, ref.full, places=5)
        self.assertAlmostEqual(result.loglik_null, ref['null'], places=5)

    def test_filter_precedes_fitting_and_fdr(self):
        from statsmodels.stats.multitest import multipletests
        stats = pd.read_csv(self.root / "out" / "component_statistics.tsv", sep="\t")
        self.assertNotIn("SE_4", stats.event_id.tolist())
        self.assertIn("SE_5", stats.event_id.tolist())  # Covered boundary response stays in the FDR family.
        filtered = self.events.status.eq("filtered_low_coverage")
        passing = self.events.loc[~filtered]
        expected = multipletests(passing.p_event.fillna(1), method="fdr_bh")[1]
        np.testing.assert_allclose(passing.q_event.dropna(), expected[passing.p_event.notna()])
        self.assertTrue(self.events.loc[filtered, ["p_event", "q_event"]].isna().all().all())
        manifest = json.loads((self.root / "out" / "analysis.json").read_text())
        self.assertEqual(manifest["n_events_passed"], 4)
        self.assertEqual(manifest["n_events_filtered"], 1)
        self.assertEqual(manifest["filter_group_columns"], ["genotype", "treatment"])

    def test_coverage_loss_of_design_cell_filters_every_event(self):
        counts = pd.read_csv(self.root / "junctions.bed", sep="\t")
        missing_cell = self.metadata[(self.metadata.genotype == "KO") & (self.metadata.treatment == "Drug")]['sample']
        counts.loc[:, missing_cell] = 0
        counts.to_csv(self.root / "missing_cell.bed", sep="\t", index=False)
        command = self.base.copy()
        command[2] = str(self.root / "missing_cell.bed")
        command[4] = str(self.root / "missing_cell")
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        stats = pd.read_csv(self.root / "missing_cell" / "component_statistics.tsv", sep="\t")
        self.assertTrue(stats.empty)
        results = pd.read_csv(self.root / "missing_cell" / "event_statistics.tsv", sep="\t")
        self.assertEqual(len(results), 5)
        self.assertTrue(results.status.eq("filtered_low_coverage").all())
        self.assertTrue(results[["p_event", "q_event"]].isna().all().all())
        audit = pd.read_csv(self.root / "missing_cell" / "event_filter.tsv", sep="\t")
        self.assertFalse(audit.filter_pass.any())
        self.assertEqual(len(audit), 20)
        self.assertNotIn("Fitting", (self.root / "missing_cell" / "glmmTMB.log").read_text())

    def test_empty_events_and_fractional_counts(self):
        empty = self.root / "empty_events"
        empty.mkdir()
        for kind in EVENT_TYPES:
            pd.read_csv(self.root / "events" / f"EVENT_{kind}.txt", sep="\t").iloc[:0].to_csv(empty / f"EVENT_{kind}.txt", sep="\t", index=False)
        command = self.base.copy(); command[3] = str(empty); command[4] = str(self.root / "empty_out")
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertTrue(pd.read_csv(self.root / "empty_out" / "event_statistics.tsv", sep="\t").empty)
        counts = pd.read_csv(self.root / "junctions.bed", sep="\t")
        counts[counts.columns[4]] = counts[counts.columns[4]].astype(float)
        counts.iloc[0, 4] = .5
        counts.to_csv(self.root / "fractional.bed", sep="\t", index=False)
        command[2] = str(self.root / "fractional.bed")
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("integers", run.stderr)

    def test_reverse_contrast_parallel_and_sample_order(self):
        path = self.root / "contrasts.json"
        path.write_text(json.dumps({"reverse": {"genotypeKO:treatmentDrug": -1},
                                    "drug_in_KO": {"treatmentDrug": 1, "genotypeKO:treatmentDrug": 1}}))
        meta = self.metadata.sample(frac=1, random_state=2)
        meta.to_csv(self.root / "shuffled.tsv", sep="\t", index=False)
        command = self.base[:-2] + ["--contrast-file", str(path), "-p", "2"]
        command[4] = str(self.root / "reverse")
        command[command.index("--sample-metadata") + 1] = str(self.root / "shuffled.tsv")
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        frame = pd.read_csv(self.root / "reverse" / "event_statistics.tsv", sep="\t")
        reverse = frame[frame.contrast == "reverse"].set_index("event_id")
        original = self.events.set_index("event_id")
        np.testing.assert_allclose(reverse.p_event, original.p_event, atol=1e-6, equal_nan=True)
        np.testing.assert_allclose(reverse.estimate_min, -original.estimate_max, atol=1e-6, equal_nan=True)
        self.assertEqual(frame[(frame.contrast == "drug_in_KO") & (frame.event_id == "SE_1")]["Diff events"].iloc[0], "Yes")

    def test_explicit_prediction_and_effect_filter(self):
        grid = self.root / "predictions.json"
        grid.write_text(json.dumps({"profiles": {
            "wc": {"genotype": "WT", "treatment": "Control"},
            "wd": {"genotype": "WT", "treatment": "Drug"},
            "kc": {"genotype": "KO", "treatment": "Control"},
            "kd": {"genotype": "KO", "treatment": "Drug"}},
            "effects": {"dd": {"wc": 1, "wd": -1, "kc": -1, "kd": 1}}}))
        command = self.base + ["--prediction-grid", str(grid), "--effect-name", "dd", "--min-effect", ".95"]
        command[4] = str(self.root / "explicit_prediction")
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        results = pd.read_csv(self.root / "explicit_prediction" / "event_statistics.tsv", sep="\t").set_index("event_id")
        original = self.events.set_index("event_id")
        np.testing.assert_allclose(results.q_event, original.q_event, equal_nan=True)
        self.assertEqual(results.loc["SE_1", "Diff events"], "No")
        predictions = pd.read_csv(self.root / "out" / "model_predictions.tsv", sep="\t")
        dd = predictions[(predictions.event_id == "SE_1") & predictions.name.str.startswith("delta_delta:")]
        self.assertAlmostEqual(results.loc["SE_1", "selected_PSI_effect"], dd.value.iloc[0])

    def test_nonestimable_design_and_unknown_coefficient(self):
        meta = self.metadata.copy()
        meta["clone"] = meta.genotype
        meta.to_csv(self.root / "confounded.tsv", sep="\t", index=False)
        command = self.base + ["--design-only"]
        command[command.index("--formula") + 1] = "genotype * treatment + clone"
        command[command.index("--sample-metadata") + 1] = str(self.root / "confounded.tsv")
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("rank deficient", result.stderr)
        command = self.base[:-1] + ["not_a_coefficient", "--design-only"]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Unknown coefficient", result.stderr)

    def test_design_only_continuous_and_multi_level(self):
        meta = self.metadata.copy()
        meta["dose"] = ["low", "medium", "high"] * (len(meta) // 3)
        meta.to_csv(self.root / "dose.tsv", sep="\t", index=False)
        command = self.base[:-2] + ["--design-only"]
        command[4] = str(self.root / "design_only")
        command[command.index("--sample-metadata") + 1] = str(self.root / "dose.tsv")
        command[command.index("--formula") + 1] = "genotype * treatment + dose + time"
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        design = pd.read_csv(self.root / "design_only" / "design_matrix.tsv", sep="\t")
        self.assertIn("time", design)
        self.assertEqual(len([c for c in design if c.startswith("dose")]), 2)

    def test_all_event_types(self):
        directory = self.root / "all_types"
        directory.mkdir()
        (directory / "events").mkdir()
        counts = pd.read_csv(self.root / "junctions.bed", sep="\t").iloc[:3].copy()
        values = counts.iloc[:, 4:].to_numpy()
        counts = []
        for index, kind in enumerate(EVENT_TYPES, 1):
            row = pd.read_csv(ROOT / "test" / "data" / f"EVENT_{kind}.txt", sep="\t").iloc[0].copy()
            row["event_id"] = f"{kind}_1"
            # Existing fixture coordinates can overlap between event types; use
            # separate chromosomes while retaining each type's geometry.
            for column in row.index:
                if isinstance(row[column], str):
                    row[column] = row[column].replace("chr10", f"chr{index + 20}")
            pd.DataFrame([row]).to_csv(directory / "events" / f"EVENT_{kind}.txt", sep="\t", index=False)
            pairs = component_pairs(kind, row)
            for side, ids in enumerate([{a for a, _ in pairs}, {b for _, b in pairs}]):
                for key in sorted(ids):
                    chrom, coordinate = key.split(":"); start, end = coordinate.split("-")
                    counts.append([chrom, start, end, key, *values[0 if side == 0 else 2]])
        pd.DataFrame(counts, columns=["chr", "start", "end", "ID", *self.metadata['sample']]).to_csv(directory / "junctions.bed", sep="\t", index=False)
        command = self.base.copy()
        command[2:5] = [str(directory / "junctions.bed"), str(directory / "events"), str(directory / "out")]
        command += ["--excel", "True", "-p", "2"]
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        result = pd.read_csv(directory / "out" / "event_statistics.tsv", sep="\t")
        self.assertEqual(set(result.event_type), set(EVENT_TYPES))
        self.assertTrue(result["Diff events"].eq("Yes").all(), result[["event_type", "status", "p_event"]])
        self.assertTrue((directory / "out" / "factorial_results.xlsx").exists())

    def test_report(self):
        from lib.factorial_report import write_report
        root = self.root / "report_input"
        root.mkdir(exist_ok=True)
        if not (root / "splicing").exists():
            (root / "splicing").symlink_to(self.root / "out")
        write_report(root, self.root / "plots")
        self.assertTrue((self.root / "plots" / "summary.html").exists())
        self.assertTrue((self.root / "plots" / "png" / "barplot_splicing_summary.png").exists())


if __name__ == "__main__":
    unittest.main()
