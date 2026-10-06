"""Event-centric beta-binomial orchestration; regression is owned by glmmTMB."""
import csv
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

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
    """Conjunction test; only coverage-passing events enter the FDR family."""
    rows = []
    passed = (None if filtering is None else
              set(filtering.loc[filtering.filter_pass.astype(bool), KEYS].itertuples(index=False, name=None)))
    grouped = {(a, b, c): g for (a, b, c), g in statistics.groupby(KEYS + ["contrast"], sort=False)}
    for (kind, event), components in definitions.groupby(KEYS, sort=False):
        expected = set(components.component_id)
        filtered = passed is not None and (kind, event) not in passed
        for contrast in contrasts:
            g = grouped.get((kind, event, contrast), pd.DataFrame(columns=STAT_COLUMNS))
            valid = (len(g) == len(expected) and set(g.component_id) == expected and
                     g.status.eq("ok").all() and np.isfinite(g.p_component).all() and
                     g.p_component.between(0, 1).all() and np.isfinite(g.estimate).all())
            same_direction = valid and (g.estimate.gt(0).all() or g.estimate.lt(0).all())
            status = "filtered_low_coverage" if filtered else "ok" if same_direction else "inconsistent_direction" if valid else "untestable"
            estimate = g.estimate.dropna()
            rows.append(dict(event_type=kind, event_id=event, contrast=contrast,
                             n_components=len(expected), n_components_ok=int(g.status.eq("ok").sum()),
                             status=status, p_event=np.nan if filtered else float(g.p_component.max()) if same_direction else 1.0 if valid else np.nan,
                             estimate_min=estimate.min(), estimate_max=estimate.max(),
                             conservative_component_effect=(float(estimate.iloc[np.argmin(abs(estimate))]) if same_direction else np.nan),
                             n_samples=int(g.n_samples.min()) if len(g) else 0))
    columns = KEYS + ["contrast", "n_components", "n_components_ok", "status", "p_event",
                      "estimate_min", "estimate_max", "conservative_component_effect", "n_samples"]
    result = pd.DataFrame(rows, columns=columns)
    result["q_event"] = np.nan
    for _, group in result.groupby("contrast", sort=False):
        group = group.loc[group.status.ne("filtered_low_coverage")]
        if group.empty:
            continue
        # Do not silently remove failed/untestable events from the correction family.
        q = multipletests(group.p_event.fillna(1), method="fdr_bh" if adjust == "BH" else "fdr_by")[1]
        result.loc[group.index, "q_event"] = np.where(group.p_event.notna(), q, np.nan)
    result["Diff events"] = np.where(result.status.eq("ok") & result.q_event.lt(fdr), "Yes", "No")
    return result


def aggregate_predictions(predictions, definitions):
    expected = definitions.groupby(KEYS).size()
    rows = []
    for key, group in predictions.groupby(KEYS + ["kind", "name"], sort=False):
        complete = (len(group) == expected.loc[key[:2]] and not group.component_id.duplicated().any()
                    and np.isfinite(group.value).all())
        rows.append([*key, group.value.mean() if complete else np.nan,
                     group.value.min() if complete else np.nan, group.value.max() if complete else np.nan,
                     "ok" if complete else "incomplete_components"])
    return pd.DataFrame(rows, columns=KEYS + ["kind", "name", "value", "component_min", "component_max", "status"])


def observed_psi(events, junctions, samples, minimum_reads):
    functions = {"SE": (shibalib.se, shibalib.col_se), "MSE": (shibalib.mse, shibalib.col_mse),
                 "MXE": (shibalib.mxe, shibalib.col_mxe), "RI": (shibalib.ri, shibalib.col_ri)}
    matrices, tables = [], {}
    for kind in EVENT_TYPES:
        function, cols = functions.get(kind, (shibalib.afe_ale if kind in ("AFE", "ALE") else shibalib.five_three,
                                             shibalib.col_five_three_afe_ale))
        table = shibalib.make_psi_table_sample(samples, events[kind], junctions, function, cols, 1, minimum_reads)
        tables[kind] = table
        matrix = table[["event_id", "pos_id"] + [s + "_PSI" for s in samples]].copy()
        matrix.columns = ["event_id", "pos_id"] + samples
        matrices.append(matrix)
    return pd.concat(matrices, ignore_index=True), tables


def run(args):
    validate_options(args)
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
    log = output / "glmmTMB.log"
    with tempfile.TemporaryDirectory(prefix=".factorial-", dir=output) as directory:
        work = Path(directory)
        stage = work / "results"; stage.mkdir()
        meta.to_csv(work / "metadata.tsv", sep="\t", index=False)
        cfg = dict(output=str(stage), metadata=str(work / "metadata.tsv"), formula=args.formula,
                   coefficients=args.coef, contrasts=contrasts, reference_levels=references,
                   categorical=args.categorical, continuous=args.continuous, prediction=prediction,
                   minimum_reads=args.minimum_reads, processes=args.num_process, design_only=args.design_only,
                   design_rds=str(work / "design.rds"), batches=[])
        config_file = work / "config.json"
        config_file.write_text(json.dumps(cfg, allow_nan=False))
        run_r(args, "design", config_file, log)
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
                if counts.empty:
                    continue
                path = work / f"{kind}_{start}.tsv"
                counts.to_csv(path, sep="\t", index=False)
                cfg["batches"].append(str(path))
        definitions = pd.concat(definitions, ignore_index=True) if definitions else pd.DataFrame(columns=COMPONENT_COLUMNS)
        filtering = pd.concat(filter_tables, ignore_index=True) if filter_tables else pd.DataFrame(columns=FILTER_COLUMNS)
        filtering.to_csv(stage / "event_filter.tsv", sep="\t", index=False)
        n_passed = len(filtering.loc[filtering.filter_pass.astype(bool), KEYS].drop_duplicates())
        n_events = len(definitions[KEYS].drop_duplicates())
        logger.info("Coverage filter: %d / %d events passed (reads >= %d; fraction >= %g in every group)",
                    n_passed, n_events, args.minimum_reads, args.min_sample_fraction)
        definitions.to_csv(stage / "components.tsv", sep="\t", index=False)
        config_file.write_text(json.dumps(cfg, allow_nan=False))
        if cfg["batches"]:
            run_r(args, "fit", config_file, log)
        statistics = (pd.read_csv(stage / "component_statistics.tsv", sep="\t", keep_default_na=False, na_values=["NA"])
                      if (stage / "component_statistics.tsv").exists() else pd.DataFrame(columns=STAT_COLUMNS))
        predictions = (pd.read_csv(stage / "component_predictions.tsv", sep="\t", keep_default_na=False, na_values=["NA"])
                       if (stage / "component_predictions.tsv").exists() else pd.DataFrame(columns=PRED_COLUMNS))
        statistics.to_csv(stage / "component_statistics.tsv", sep="\t", index=False)
        predictions.to_csv(stage / "component_predictions.tsv", sep="\t", index=False)
        tests = args.coef + list(contrasts)
        result = aggregate_events(statistics, definitions, tests, args.fdr, args.p_adjust, filtering)
        fitted = aggregate_predictions(predictions, definitions)
        fitted.to_csv(stage / "model_predictions.tsv", sep="\t", index=False)
        if args.min_effect:
            effect = fitted[(fitted.kind == "effect") & (fitted.name == args.effect_name)][KEYS + ["value"]]
            result = result.merge(effect.rename(columns={"value": "selected_PSI_effect"}), on=KEYS, how="left")
            result.loc[~result.selected_PSI_effect.abs().ge(args.min_effect), "Diff events"] = "No"
        matrix, observed = observed_psi(events, junctions, samples, args.minimum_reads)
        matrix.to_csv(stage / "PSI_matrix_sample.txt", sep="\t", index=False)
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
                    "minimum_reads": args.minimum_reads, "effect_name": args.effect_name, "min_effect": args.min_effect,
                    "event_prediction": "unweighted component mean; distinct from observed Shiba PSI",
                    "n_events": n_events, "n_samples": len(samples)}
        (stage / "analysis.json").write_text(json.dumps(manifest, indent=2))
        meta.to_csv(stage / "sample_metadata.tsv", sep="\t", index=False)
        for file in stage.iterdir():
            os.replace(file, output / file.name)
    logger.info("Regression completed: %d event/contrast rows; %d significant; %d untestable",
                len(result), result["Diff events"].eq("Yes").sum(), result.status.eq("untestable").sum())
