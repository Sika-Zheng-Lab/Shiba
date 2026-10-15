"""Event-centric beta-binomial orchestration; regression is owned by glmmTMB."""
import csv
from concurrent.futures import ProcessPoolExecutor
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from . import shibalib
from .factorial_config import config_arguments
from .event_components import EVENT_TYPES, COMPONENT_COLUMNS, extract_components

logger = logging.getLogger(__name__)
KEYS = ["event_type", "event_id"]
STAT_COLUMNS = KEYS + ["component_id", "contrast", "estimate", "se", "ci_low", "ci_high",
                       "p_component", "lr", "rho", "phi", "loglik_full", "loglik_null",
                       "n_samples", "full_convergence", "full_pd_hessian", "null_convergence",
                       "null_pd_hessian", "status", "diagnostic"]
PRED_COLUMNS = KEYS + ["component_id", "kind", "name", "value"]
FILTER_COLUMNS = KEYS + ["filter_group", "n_samples", "n_eligible", "n_required", "group_pass", "filter_pass"]


class StageTimer:
    def __init__(self):
        self.last = time.perf_counter()
        self.rows = []

    def mark(self, name):
        now = time.perf_counter()
        self.rows.append({"stage": name, "elapsed_seconds": now - self.last})
        self.last = now
        logger.info("Stage %s: %.3f seconds", name, self.rows[-1]["elapsed_seconds"])


class PassingBatchWriter:
    """Bounded batches of passing events, independent of input batch/type boundaries."""
    def __init__(self, directory, event_limit):
        self.directory, self.event_limit = Path(directory), event_limit
        self.frames, self.paths, self.n_events = [], [], 0

    def add(self, counts):
        if counts.empty:
            return
        numbers = counts.groupby(KEYS, sort=False).ngroup().to_numpy()
        total = int(numbers.max()) + 1
        start = 0
        while start < total:
            stop = min(total, start + self.event_limit - self.n_events)
            self.frames.append(counts.loc[(numbers >= start) & (numbers < stop)])
            self.n_events += stop - start
            start = stop
            if self.n_events == self.event_limit:
                self.flush()

    def flush(self):
        if self.frames:
            path = self.directory / f"passing_{len(self.paths):06d}.tsv"
            pd.concat(self.frames, ignore_index=True).to_csv(path, sep="\t", index=False)
            self.paths.append(str(path))
            self.frames, self.n_events = [], 0


def add_arguments(parser):
    parser.add_argument("--stat-method", choices=["legacy", "beta-binomial"], default="legacy")
    parser.add_argument("--sample-metadata", help="TSV with sample and formula variables; defaults to --group")
    parser.add_argument("--formula", help="R fixed-effect formula, e.g. genotype * treatment + batch")
    parser.add_argument("--coef", action="append", default=[], help="R design coefficient name; repeat for multiple tests")
    parser.add_argument("--contrast-file", help="JSON mapping contrast names to {coefficient: weight}")
    parser.add_argument("--reference-level", action="append", default=[], metavar="FACTOR=LEVEL")
    parser.add_argument("--categorical", action="append", default=[], help="Force metadata variable to be categorical")
    parser.add_argument("--continuous", action="append", default=[], help="Force metadata variable to be numeric")
    parser.add_argument("--prediction-grid", help="JSON with named profiles and PSI effect weights")
    parser.add_argument("--effect-name", help="PSI effect used with --min-effect")
    parser.add_argument("--min-effect", type=float, default=0.0, help="Optional absolute PSI effect filter, after FDR")
    parser.add_argument("--p-adjust", choices=["BH", "BY"], default="BH")
    parser.add_argument("--design-only", action="store_true", help="Validate/export design and coefficients without fitting")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable with glmmTMB and jsonlite")
    parser.add_argument("--event-batch-size", type=int, default=200)
    parser.add_argument("--component-cache-size", type=int, default=10000,
                        help="Maximum exact component results cached per run; 0 disables caching")
    parser.add_argument("--min-sample-fraction", type=float, default=0.5,
                        help="Required fraction of covered samples in every filter group (default: 0.5)")
    parser.add_argument("--filter-group", action="append", default=[], metavar="COLUMN",
                        help="Metadata column defining coverage groups; repeat for combinations. Default: categorical design factors")


def validate_options(args):
    if not args.formula or not (args.sample_metadata or args.group):
        raise ValueError("beta-binomial requires --formula and --sample-metadata (or --group)")
    if args.onlypsi or args.onlypsi_group or args.beta_binomial or args.ttest:
        raise ValueError("beta-binomial regression cannot be combined with legacy -b/-t or PSI-only modes")
    if not args.design_only and not (args.coef or args.contrast_file):
        raise ValueError("Specify --coef or --contrast-file, or use --design-only to list coefficients")
    if args.num_process < 1 or args.event_batch_size < 1 or args.minimum_reads < 1:
        raise ValueError("Processes, event batch size and minimum reads must be positive")
    if args.component_cache_size < 0:
        raise ValueError("--component-cache-size must be nonnegative")
    if not 0 < args.fdr < 1 or not np.isfinite(args.min_effect) or args.min_effect < 0:
        raise ValueError("Require 0 < FDR < 1 and a nonnegative finite minimum effect")
    if args.min_effect and not args.effect_name:
        raise ValueError("--min-effect requires --effect-name")
    if not np.isfinite(args.min_sample_fraction) or not 0 < args.min_sample_fraction <= 1:
        raise ValueError("--min-sample-fraction must be finite and in (0, 1]")
    if len(set(args.coef)) != len(args.coef):
        raise ValueError("Duplicate coefficient tests")


def read_json(path):
    if not path:
        return None
    def unique_object(pairs):
        obj = {}
        for key, value in pairs:
            if key in obj:
                raise ValueError(f"Duplicate JSON key: {key}")
            obj[key] = value
        return obj
    with open(path) as handle:
        return json.load(handle, object_pairs_hook=unique_object)


def read_metadata(path, samples):
    validate_header(path)
    meta = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "sample" not in meta or meta["sample"].eq("").any() or meta["sample"].duplicated().any():
        raise ValueError("Metadata requires unique nonempty sample IDs")
    if set(meta["sample"]) != set(samples):
        raise ValueError("Metadata sample IDs must exactly match junction samples; subset inputs explicitly")
    return meta.set_index("sample").loc[samples].reset_index()


def validate_header(path):
    with open(path, newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    if not header or any(not column for column in header) or len(set(header)) != len(header):
        raise ValueError(f"Empty or duplicate column names in {path}")


def run_r(args, stage, config_file, log):
    executable = shutil.which(args.rscript)
    if executable is None:
        raise ValueError(f"Rscript executable not found: {args.rscript}")
    script = Path(__file__).resolve().parents[1] / "beta_binomial_glm.R"
    with open(log, "a") as handle:
        result = subprocess.run([executable, "--vanilla", str(script), stage, str(config_file)],
                                stdout=handle, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        with open(log) as handle:
            detail = handle.read()[-6000:]
        raise RuntimeError(f"glmmTMB {stage} failed (see {log}):\n{detail}")


def resolve_filter_groups(meta, design, columns):
    """Use explicit metadata columns or observed categorical design cells."""
    columns = list(columns) if columns else list(design.get("levels") or {})
    if len(set(columns)) != len(columns) or any(c not in meta or c == "sample" for c in columns):
        raise ValueError("--filter-group requires unique metadata columns other than sample")
    if columns and (meta[columns].isna().any().any() or meta[columns].eq("").any().any()):
        raise ValueError("Missing metadata in coverage filter groups")
    groups = {}
    for row in meta.to_dict("records"):
        label = json.dumps({c: row[c] for c in columns}, ensure_ascii=False, sort_keys=True)
        groups.setdefault(label, []).append(row["sample"])
    return columns, groups


def filter_events(counts, groups, minimum_reads=10, fraction=0.5):
    """Screen on the common component coverage mask, before any model fitting."""
    if counts.empty:
        return pd.DataFrame(columns=FILTER_COLUMNS)
    totals = counts.assign(total=counts.inclusion + counts.exclusion)
    coverage = totals.groupby(KEYS + ["sample"], sort=False).total.min().unstack("sample")
    rows = []
    for label, samples in groups.items():
        eligible = coverage.reindex(columns=samples).ge(minimum_reads).sum(axis=1)
        required = int(np.ceil(len(samples) * fraction))
        for (kind, event), number in eligible.items():
            rows.append([kind, event, label, len(samples), int(number), required, number >= required])
    result = pd.DataFrame(rows, columns=FILTER_COLUMNS[:-1])
    result["filter_pass"] = result.groupby(KEYS, sort=False).group_pass.transform("all")
    return result


def aggregate_events(statistics, definitions, contrasts, fdr=.05, adjust="BH", filtering=None):
    """Vectorized conjunction tests; retain post-filter failures in each FDR family."""
    ids = KEYS + ["contrast"]
    expected = definitions.groupby(KEYS, sort=False).component_id.nunique().rename("n_components").reset_index()
    result = expected.merge(pd.DataFrame({"contrast": contrasts}), how="cross")
    stats = statistics[ids + ["component_id", "status", "estimate", "p_component", "n_samples"]].merge(
                             definitions[KEYS + ["component_id"]].drop_duplicates().assign(known=True),
                             on=KEYS + ["component_id"], how="left", validate="many_to_one")
    for column in ["estimate", "p_component", "n_samples"]:
        stats[column] = pd.to_numeric(stats[column])
    stats["ok"] = stats.status.eq("ok")
    stats["valid"] = (stats.ok & stats.known.eq(True) & np.isfinite(stats.estimate) &
                      np.isfinite(stats.p_component) & stats.p_component.between(0, 1))
    metrics = stats.groupby(ids, sort=False).agg(
        rows=("component_id", "size"), distinct=("component_id", "nunique"),
        valid=("valid", "all"), n_components_ok=("ok", "sum"),
        p_max=("p_component", "max"), estimate_min=("estimate", "min"),
        estimate_max=("estimate", "max"), n_samples=("n_samples", "min")).reset_index()
    result = result.merge(metrics, on=ids, how="left", validate="one_to_one")
    valid = (result.valid.eq(True) & result.rows.eq(result.n_components) &
             result.distinct.eq(result.n_components))
    direction = valid & (result.estimate_min.gt(0) | result.estimate_max.lt(0))
    filtered = pd.Series(False, index=result.index)
    if filtering is not None:
        passed = filtering.loc[filtering.filter_pass.astype(bool), KEYS].drop_duplicates().assign(passed=True)
        result = result.merge(passed, on=KEYS, how="left", validate="many_to_one")
        filtered = ~result.passed.eq(True)
    result["status"] = np.select([filtered, direction, valid],
                                  ["filtered_low_coverage", "ok", "inconsistent_direction"], default="untestable")
    result["p_event"] = np.where(filtered, np.nan, np.where(direction, result.p_max, np.where(valid, 1., np.nan)))
    result["conservative_component_effect"] = np.where(direction,
        np.where(result.estimate_min.gt(0), result.estimate_min, result.estimate_max), np.nan)
    for column in ["n_components_ok", "n_samples"]:
        result[column] = result[column].fillna(0).astype(int)
    result["q_event"] = np.nan
    for _, group in result.loc[~filtered].groupby("contrast", sort=False):
        q = multipletests(group.p_event.fillna(1), method="fdr_bh" if adjust == "BH" else "fdr_by")[1]
        result.loc[group.index, "q_event"] = np.where(group.p_event.notna(), q, np.nan)
    result["Diff events"] = np.where(result.status.eq("ok") & result.q_event.lt(fdr), "Yes", "No")
    return result[KEYS + ["contrast", "n_components", "n_components_ok", "status", "p_event",
                         "estimate_min", "estimate_max", "conservative_component_effect", "n_samples", "q_event", "Diff events"]]


def aggregate_predictions(predictions, definitions):
    ids = KEYS + ["kind", "name"]
    columns = ids + ["value", "component_min", "component_max", "status"]
    if predictions.empty:
        return pd.DataFrame(columns=columns)
    expected = definitions.groupby(KEYS).size().rename("expected").reset_index()
    result = predictions.assign(finite=np.isfinite(predictions.value)).groupby(ids, sort=False).agg(
        rows=("component_id", "size"), distinct=("component_id", "nunique"),
        finite=("finite", "all"), value=("value", "mean"),
        component_min=("value", "min"), component_max=("value", "max")).reset_index()
    result = result.merge(expected, on=KEYS, how="left", validate="many_to_one")
    complete = result.rows.eq(result.expected) & result.distinct.eq(result.rows) & result.finite
    result.loc[~complete, ["value", "component_min", "component_max"]] = np.nan
    result["status"] = np.where(complete, "ok", "incomplete_components")
    return result[columns]


def observed_psi(events, junctions, samples, minimum_reads, num_process=1):
    functions = {"SE": (shibalib.se, shibalib.col_se), "MSE": (shibalib.mse, shibalib.col_mse),
                 "MXE": (shibalib.mxe, shibalib.col_mxe), "RI": (shibalib.ri, shibalib.col_ri)}
    matrices, tables = [], {}
    shm, shm_info = None, None
    executor = None
    try:
        if num_process > 1 and isinstance(junctions, shibalib.JunctionData) and junctions.array.nbytes:
            shm, shm_info = junctions.to_shared_memory()
            workers = min(num_process, max(len(frame) for frame in events.values()))
            if workers > 1:
                executor = ProcessPoolExecutor(
                    max_workers=workers, initializer=shibalib._init_junc_worker,
                    initargs=(shm_info["shm_name"], shm_info["shape"],
                              shm_info["junction_ids"], shm_info["sample_ids"]))
        for kind in EVENT_TYPES:
            function, cols = functions.get(kind, (shibalib.afe_ale if kind in ("AFE", "ALE") else shibalib.five_three,
                                                 shibalib.col_five_three_afe_ale))
            table = shibalib.make_psi_table_sample(samples, events[kind], junctions, function, cols, num_process, minimum_reads, shm_info=shm_info, executor=executor)
            tables[kind] = table
            matrix = table[["event_id", "pos_id"] + [s + "_PSI" for s in samples]].copy()
            matrix.columns = ["event_id", "pos_id"] + samples
            matrices.append(matrix)
    finally:
        try:
            if executor is not None:
                executor.shutdown(wait=True)
        finally:
            if shm is not None:
                shm.close()
                shm.unlink()
    return pd.concat(matrices, ignore_index=True), tables


def run(args):
    validate_options(args)
    timer = StageTimer()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    logger.info("Running event-centric beta-binomial regression using glmmTMB")
    # Check before JunctionData's integer conversion, which would truncate fractions.
    validate_header(args.junctions)
    raw = pd.read_csv(args.junctions, sep="\t", dtype=str, keep_default_na=False)
    if list(raw.columns[:4]) != ["chr", "start", "end", "ID"] or raw.ID.duplicated().any():
        raise ValueError("Junction input requires chr/start/end/ID and unique junction IDs")
    samples = list(raw.columns[4:])
    if not samples:
        raise ValueError("No junction samples")
    values = raw.iloc[:, 4:].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    if np.any(~np.isfinite(values) | (values < 0) | (values != np.floor(values)) | (values >= 2**53)):
        raise ValueError("Junction counts must be nonnegative, exactly representable integers below 2^53")
    junctions = shibalib.JunctionData.from_dataframe(raw)
    del raw, values
    meta = read_metadata(args.sample_metadata or args.group, samples)
    references = {}
    for entry in args.reference_level:
        factor, separator, level = entry.partition("=")
        if not separator or not factor or not level or factor in references:
            raise ValueError("Reference levels must be unique FACTOR=LEVEL entries")
        references[factor] = level
    contrasts = read_json(args.contrast_file) or {}
    if not isinstance(contrasts, dict) or any(not isinstance(v, dict) for v in contrasts.values()):
        raise ValueError("Contrast file must map test names to coefficient-weight objects")
    prediction = read_json(args.prediction_grid)
    timer.mark("input_validation")
    log = output / "glmmTMB.log"
    with tempfile.TemporaryDirectory(prefix=".factorial-", dir=output) as directory:
        work = Path(directory)
        stage = work / "results"; stage.mkdir()
        meta.to_csv(work / "metadata.tsv", sep="\t", index=False)
        cfg = dict(output=str(stage), metadata=str(work / "metadata.tsv"), formula=args.formula,
                   coefficients=args.coef, contrasts=contrasts, reference_levels=references,
                   categorical=args.categorical, continuous=args.continuous, prediction=prediction,
                   minimum_reads=args.minimum_reads, processes=args.num_process, design_only=args.design_only,
                   design_rds=str(work / "design.rds"), batches=[], component_cache_size=args.component_cache_size)
        config_file = work / "config.json"
        config_file.write_text(json.dumps(cfg, allow_nan=False))
        run_r(args, "design", config_file, log)
        timer.mark("design")
        design = read_json(stage / "design.json")
        filter_columns, filter_groups = resolve_filter_groups(meta, design, args.filter_group)
        effect_names = ((design.get("prediction") or {}).get("effects") or {}).keys()
        if args.effect_name and args.effect_name not in effect_names:
            raise ValueError(f"Unknown PSI effect {args.effect_name}; available: {list(effect_names)}")
        if args.design_only:
            if (output / "analysis.json").exists():
                raise ValueError("Use a separate output directory for --design-only; this directory contains an analysis")
            for file in stage.iterdir():
                os.replace(file, output / file.name)
            logger.info("Design exported; available coefficients: %s", design["coefficients"])
            return
        events = shibalib.read_events(args.event)
        definitions, filter_tables = [], []
        batch_writer = PassingBatchWriter(work, max(args.event_batch_size, 4 * args.num_process))
        for kind in EVENT_TYPES:
            if events[kind].event_id.isna().any() or events[kind].event_id.duplicated().any():
                raise ValueError(f"Invalid event IDs for {kind}")
            for start in range(0, len(events[kind]), args.event_batch_size):
                batch = events[kind].iloc[start:start + args.event_batch_size]
                counts, component = extract_components(kind, batch, junctions, samples)
                filtering = filter_events(counts, filter_groups, args.minimum_reads, args.min_sample_fraction)
                filter_tables.append(filtering)
                passed_ids = filtering.loc[filtering.filter_pass, "event_id"].unique()
                counts = counts.loc[counts.event_id.isin(passed_ids)]
                definitions.append(component)
                batch_writer.add(counts)
        batch_writer.flush()
        cfg["batches"] = batch_writer.paths
        definitions = pd.concat(definitions, ignore_index=True) if definitions else pd.DataFrame(columns=COMPONENT_COLUMNS)
        filtering = pd.concat(filter_tables, ignore_index=True) if filter_tables else pd.DataFrame(columns=FILTER_COLUMNS)
        filtering.to_csv(stage / "event_filter.tsv", sep="\t", index=False)
        n_passed = len(filtering.loc[filtering.filter_pass.astype(bool), KEYS].drop_duplicates())
        n_events = len(definitions[KEYS].drop_duplicates())
        logger.info("Coverage filter: %d / %d events passed (reads >= %d; fraction >= %g in every group)",
                    n_passed, n_events, args.minimum_reads, args.min_sample_fraction)
        definitions.to_csv(stage / "components.tsv", sep="\t", index=False)
        config_file.write_text(json.dumps(cfg, allow_nan=False))
        timer.mark("coverage_and_batching")
        if cfg["batches"]:
            run_r(args, "fit", config_file, log)
        else:
            (stage / "backend_performance.json").write_text(json.dumps({
                "processes": args.num_process, "scheduling": "no passing events",
                "component_cache_size": args.component_cache_size, "components": 0,
                "cache_hits": 0, "optimizer_attempts": 0, "retries": 0, "elapsed_seconds": 0}, indent=2))
            pd.DataFrame(columns=KEYS + ["component_id", "cache_hit", "worker_pid", "elapsed_seconds",
                         "full_seconds", "null_seconds", "full_attempts", "null_attempts", "retries"]).to_csv(
                             stage / "component_timings.tsv", sep="\t", index=False)
            pd.DataFrame(columns=["batch", "components", "fitted_components", "cache_hits", "full_attempts",
                                  "null_attempts", "retries", "elapsed_seconds"]).to_csv(
                                      stage / "batch_timings.tsv", sep="\t", index=False)
        timer.mark("backend_fit")
        statistics = (pd.read_csv(stage / "component_statistics.tsv", sep="\t", keep_default_na=False, na_values=["NA"])
                      if (stage / "component_statistics.tsv").exists() else pd.DataFrame(columns=STAT_COLUMNS))
        predictions = (pd.read_csv(stage / "component_predictions.tsv", sep="\t", keep_default_na=False, na_values=["NA"])
                       if (stage / "component_predictions.tsv").exists() else pd.DataFrame(columns=PRED_COLUMNS))
        if not (stage / "component_statistics.tsv").exists():
            statistics.to_csv(stage / "component_statistics.tsv", sep="\t", index=False)
        if not (stage / "component_predictions.tsv").exists():
            predictions.to_csv(stage / "component_predictions.tsv", sep="\t", index=False)
        timer.mark("read_component_results")
        tests = args.coef + list(contrasts)
        result = aggregate_events(statistics, definitions, tests, args.fdr, args.p_adjust, filtering)
        timer.mark("aggregate_events")
        fitted = aggregate_predictions(predictions, definitions)
        fitted.to_csv(stage / "model_predictions.tsv", sep="\t", index=False)
        timer.mark("aggregate_predictions_and_write")
        if args.min_effect:
            effect = fitted[(fitted.kind == "effect") & (fitted.name == args.effect_name)][KEYS + ["value"]]
            result = result.merge(effect.rename(columns={"value": "selected_PSI_effect"}), on=KEYS, how="left")
            result.loc[~result.selected_PSI_effect.abs().ge(args.min_effect), "Diff events"] = "No"
        matrix, observed = observed_psi(events, junctions, samples, args.minimum_reads, args.num_process)
        matrix.to_csv(stage / "PSI_matrix_sample.txt", sep="\t", index=False)
        timer.mark("observed_psi_and_matrix")
        annotated = []
        for kind in EVENT_TYPES:
            table = result[result.event_type == kind].merge(events[kind], on="event_id", how="left", validate="many_to_one")
            if args.individual_psi:
                table = table.merge(observed[kind][["event_id"] + [s + "_PSI" for s in samples]], on="event_id", how="left")
            table.to_csv(stage / f"PSI_{kind}.txt", sep="\t", index=False)
            observed[kind].to_csv(stage / f"observed_PSI_{kind}.tsv", sep="\t", index=False)
            annotated.append(table)
        result = pd.concat(annotated, ignore_index=True)
        result.to_csv(stage / "event_statistics.tsv", sep="\t", index=False)
        summary = result.groupby(["event_type", "contrast", "status", "Diff events"], dropna=False).size().reset_index(name="Number")
        summary.to_csv(stage / "summary.txt", sep="\t", index=False)
        if args.excel:
            with pd.ExcelWriter(stage / "factorial_results.xlsx") as writer:
                for kind, table in zip(EVENT_TYPES, annotated):
                    table.to_excel(writer, sheet_name=kind, index=False)
                fitted.to_excel(writer, sheet_name="model_predictions", index=False)
        manifest = {"stat_method": "beta-binomial", "backend": "glmmTMB", "formula": args.formula,
                    "contrasts": tests, "p_adjust": args.p_adjust, "fdr": args.fdr,
                    "fdr_family": "coverage-passing events across event types, separately per contrast; post-filter failures counted as p=1",
                    "min_sample_fraction": args.min_sample_fraction, "filter_group_columns": filter_columns,
                    "filter_groups": filter_groups, "n_events_passed": n_passed, "n_events_filtered": n_events - n_passed,
                    "processes": args.num_process, "component_cache_size": args.component_cache_size,
                    "passing_batch_event_limit": batch_writer.event_limit,
                    "minimum_reads": args.minimum_reads, "effect_name": args.effect_name, "min_effect": args.min_effect,
                    "event_prediction": "unweighted component mean; distinct from observed Shiba PSI",
                    "n_events": n_events, "n_samples": len(samples)}
        (stage / "analysis.json").write_text(json.dumps(manifest, indent=2))
        meta.to_csv(stage / "sample_metadata.tsv", sep="\t", index=False)
        timer.mark("write_event_results")
        pd.DataFrame(timer.rows).to_csv(stage / "stage_timings.tsv", sep="\t", index=False)
        for file in stage.iterdir():
            os.replace(file, output / file.name)
    logger.info("Regression completed: %d event/contrast rows; %d significant; %d untestable",
                len(result), result["Diff events"].eq("Yes").sum(), result.status.eq("untestable").sum())
