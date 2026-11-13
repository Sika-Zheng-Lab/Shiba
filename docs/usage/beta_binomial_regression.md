# Event-centric beta-binomial regression

Shiba and SnakeShiba support an optional `stat_method: beta_binomial` mode for
factorial differential splicing. Python extracts integer junction comparisons;
R/glmmTMB fits each component across biological replicates. The default `fisher`
mode runs Fisher testing, with an optional Welch t-test. The old supplementary
`beta_binomial` option and `psi.py -b/--beta-binomial` have been removed.

## Installation

Use the checked-out Shiba source containing this backend, with R packages
`glmmTMB` and `jsonlite` installed in the R environment used by `Rscript`.
For an existing Conda environment:

```bash
conda activate shiba
conda install -c conda-forge r-glmmtmb r-jsonlite
Rscript --vanilla -e 'stopifnot(requireNamespace("glmmTMB", quietly=TRUE), requireNamespace("jsonlite", quietly=TRUE))'
```

See the [installation guide](../installation.md) for new Shiba/MameShiba
environments, CRAN and user-library installation, verification, and Docker or
SnakeShiba container setup. Installing the packages alone does not add this
mode to an older Shiba executable; run `python ./src/psi.py` or
`python ./shiba.py` from the matching source checkout.

`--rscript /path/to/Rscript` selects another R installation; the corresponding
pipeline setting is `rscript`. The backend uses `Rscript --vanilla`, so a custom
library must be visible through standard R library paths or `R_LIBS_USER`.
Shiba records R and glmmTMB versions in `design.json`. No Patsy or Formulaic is
used by the new engine.

Matching [sample metadata](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/beta_binomial/samples.tsv),
[experiment table](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/beta_binomial/Shiba_experiment.tsv),
[contrasts](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/beta_binomial/contrasts.json), and
[prediction profiles](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/beta_binomial/predictions.json) are provided
for the dedicated config. See the [factorial quick start](../quickstart/beta_binomial.md)
for copying, editing, and running them. These templates require your own BAMs
and reference GTF.

## Inputs and model

Use the existing `junctions.bed` and `EVENT_*.txt` files. Supply metadata with one
row per independent biological replicate and exactly the same sample IDs as the
junction table. Row order does not matter. Missing covariates, duplicate sample
IDs and rank-deficient designs are errors; samples are not silently discarded.

```text
sample  genotype  treatment  batch
WT_C_1  WT        Control    B1
WT_C_2  WT        Control    B2
WT_D_1  WT        Drug       B1
WT_D_2  WT        Drug       B2
KO_C_1  KO        Control    B1
KO_C_2  KO        Control    B2
KO_D_1  KO        Drug       B1
KO_D_2  KO        Drug       B2
```

The actual file must be tab-separated. In full pipelines, keep the usual
`sample`, `bam`, `group` experiment table and supply a separate metadata file.
`group` remains a label for the existing pipeline, not a restriction to two
regression groups. The standalone `psi.py` regression does not require `group`.

For component c of event e and sample s:

```
K_ecs ~ BetaBinomial(N_ecs, mu_ecs, rho_ec)
N_ecs = inclusion_ecs + exclusion_ecs
logit(mu_ecs) = X_s beta_ec
```

Dispersion is estimated separately for each event/component, constant across
samples (`dispformula = ~1`). Reported `phi` is glmmTMB's concentration parameter;
`rho = 1 / (1 + phi)` is the overdispersion parameter. There is no across-event
shrinkage. This initial implementation accepts fixed-effect formulas only;
random effects, offsets, smooths and arbitrary R function calls are rejected.

Supported formula algebra includes `+`, `-`, `*`, `:`, `/`, `^`, parentheses,
and intercept removal (`~ 0 + condition`). Examples:

```
~ genotype * treatment + batch
~ condition + sex + batch
~ genotype * time
~ genotype * treatment * sex
```

Numeric columns are continuous by default; other columns become factors.
Use repeatable `--categorical` / `--continuous` options to override inference.
Factor levels are sorted, and `--reference-level FACTOR=LEVEL` sets a baseline.
Do not interpret numeric time as a categorical time course without declaring it.

## CLI

First inspect the actual R coefficient names:

```bash
python src/psi.py junctions.bed events design_check \
  --stat-method beta_binomial --sample-metadata samples.tsv \
  --formula 'genotype * treatment + batch' \
  --reference-level genotype=WT --reference-level treatment=Control \
  --design-only
```

Then test the interaction:

```bash
python src/psi.py junctions.bed events results \
  --stat-method beta_binomial --sample-metadata samples.tsv \
  --formula 'genotype * treatment + batch' \
  --reference-level genotype=WT --reference-level treatment=Control \
  --coef 'genotypeKO:treatmentDrug' -p 4
```

R coefficient names (for example `genotypeKO:treatmentDrug`) are used, rather
than Patsy-style names. Repeat `--coef` to test several coefficients. Use a
separate output directory for design inspection so that it does not replace
the design metadata associated with an existing analysis.

Arbitrary scalar linear contrasts are supplied as a JSON file:

```json
{
  "interaction": {"genotypeKO:treatmentDrug": 1},
  "drug_in_KO": {"treatmentDrug": 1, "genotypeKO:treatmentDrug": 1},
  "reversed_interaction": {"genotypeKO:treatmentDrug": -1}
}
```

Pass `--contrast-file contrasts.json`. Each contrast tests `c' beta = 0` using a
full versus constrained-null likelihood-ratio test with one degree of freedom.
The null design uses a basis for the null space of `c'`, so a general contrast
is not approximated by removing a single formula term. Multi-degree-of-freedom
omnibus tests are not part of this interface; factors with three or more levels
are supported through their coefficients and scalar contrasts.

## Component definitions and event calls

| Event | Inclusion/exclusion comparisons |
| --- | --- |
| SE | a/c and b/c |
| MSE | every inclusion junction / the final skipping junction |
| FIVE, THREE | a/b |
| MXE | a1/b1, a1/b2, a2/b1, a2/b2 |
| RI | start boundary / spliced junction; end boundary / spliced junction |
| AFE, ALE | all A-side / B-side junction pairs |

RI's positive direction means increased retention. MXE uses the current
strand-aware exon A orientation. Absent junction IDs have zero counts. Negative,
fractional, non-finite or inexactly representable input counts are rejected.
No inclusion averages or duplicated skipping counts enter the likelihood.

For each event, a sample is eligible only if **every** required component has
`inclusion + exclusion >= --minimum-reads` (default 10). All components use this
same sample set. Before fitting, the event must have eligible samples in **at
least half of every group** (`--min-sample-fraction 0.5`, the default). The
required number is rounded up: two eligible samples in a group of three, for
example. Counts are evaluated per component and sample, not pooled across
samples or junctions.

By default, groups are the **observed combinations of all categorical factors
in the model formula**, regardless of whether metadata has a `group` column.
For `genotype * treatment`, these are genotype × treatment cells. A categorical
`batch` in the formula also splits the cells; continuous covariates do not.
Numeric-coded categories must be declared using `--categorical` (or reference
levels) as usual. With no categorical factors, all samples form one group.
Override the grouping columns explicitly when needed, for example:

```bash
--filter-group genotype --filter-group treatment --min-sample-fraction 0.5
```

Only coverage-passing events are sent to glmmTMB and included in FDR correction.
Filtered events remain in event output with `status = filtered_low_coverage`,
missing P/Q values and `Diff events = No`; `event_filter.tsv` records eligible
and required sample counts per group. Observed PSI tables still contain all
input events. Coverage screening does not use observed effect direction or
P values. After screening, the remaining design must retain full rank and have more samples than
mean coefficients plus one dispersion parameter. This is an identifiability
check, not a guarantee of reliable small-sample inference.

An event requires all component fits and tests to succeed, and all tested
contrast estimates to have the same nonzero sign. Its P value is the maximum
component P value. Contradictory directions give `p_event = 1`. A missing,
boundary-only or failed component gives an untestable event; other components
cannot rescue it. Components sharing junctions are not assumed independent,
and their likelihoods are not added together.

By default BH correction pools all eight event types, separately for each
contrast. `--p-adjust BY` requests the more conservative arbitrary-dependence
correction. Only events passing the coverage filter enter each correction
family. If a passing event is subsequently untestable or fails fitting, it
contributes P=1 internally but displays missing P/Q values. Separate
contrast families do not provide FDR control over a subsequently selected union
of all contrasts. Shared junctions between events can also induce dependence;
BH's usual dependence assumptions apply.

`Diff events` uses event Q < `--fdr` and the consistency rule. Legacy odds-ratio
and `-d/--psi` cutoffs do not apply. An optional biological relevance filter uses
`--effect-name NAME --min-effect VALUE` after correction; the Q values still
refer to logit-scale tests, not an effect-size-threshold null hypothesis.

## PSI predictions and effect sizes

For a single selected contrast involving one or two binary categorical factors,
Shiba automatically generates their condition profiles. Other covariates retain
the same empirical distribution in every profile: predictions are averaged
equally over **all design samples**, not pooled by read depth. For two factors,
the first varies between strata and the second defines within-stratum changes.
The resulting PSI effects include the two changes and their difference.

For other designs, or multiple contrasts, explicitly provide
`--prediction-grid predictions.json`:

```json
{
  "profiles": {
    "WT_Control": {"genotype": "WT", "treatment": "Control"},
    "WT_Drug": {"genotype": "WT", "treatment": "Drug"},
    "KO_Control": {"genotype": "KO", "treatment": "Control"},
    "KO_Drug": {"genotype": "KO", "treatment": "Drug"}
  },
  "effects": {
    "delta_WT": {"WT_Drug": 1, "WT_Control": -1},
    "delta_KO": {"KO_Drug": 1, "KO_Control": -1},
    "delta_delta": {"KO_Drug": 1, "KO_Control": -1, "WT_Drug": -1, "WT_Control": 1}
  }
}
```

Numeric settings also support chosen time/dose values. Profile values must be
valid for the design; users are responsible for avoiding extrapolation outside
supported covariate combinations. PSI effect weights must sum to zero.

Component predictions are fitted inclusion probabilities. Event predictions
and PSI effects are the **equal-weight mean across all required components**;
they are not identical to historical Shiba PSI based on mean junction counts.
Historical observed PSI is retained separately, including its existing coverage
rule. The current implementation reports point predictions without PSI-scale
confidence intervals. Coefficient intervals are approximate Wald intervals;
P values come from the LRT.

A logit interaction of zero does not imply a PSI-scale difference-in-differences
of zero. The regression P value must not be presented as testing delta-delta PSI.
Likewise, component direction consistency is enforced on the tested logit
contrast; it does not assert equality of component effect magnitudes.

## Output

### Locations and execution modes

For `python shiba.py config.yaml` and SnakeShiba, all paths below are relative
to the configured `workdir`:

```text
workdir/
├── annotation/                  # Optional transcript assembly
├── events/                      # EVENT_<TYPE>.txt definitions
├── junctions/                   # junctions.bed and counting intermediates
├── results/
│   ├── splicing/                # All regression files listed below
│   ├── expression/              # Counts, abundance, optional DESeq2 results
│   └── pca/                     # TPM and observed-PSI PCA
├── plots/
│   ├── summary.html             # Main self-contained interactive report
│   ├── figure_style.json
│   ├── png/
│   ├── pdf/
│   ├── svg/
│   ├── interactive/             # Embedded plot fragments
│   └── data/                    # Standalone PCA HTMLs, when PCA is included
└── report.json                  # Pipeline execution provenance
```

`<TYPE>` is one of `SE`, `FIVE`, `THREE`, `MXE`, `RI`, `MSE`,
`AFE`, or `ALE`. A standalone `src/psi.py` run writes the regression files
**directly into its output argument**, without creating expression, PCA, or HTML
reports. To use the reporting command below, choose `results/splicing` as that
output and pass its parent `results` to `plots.py`.

`shiba.py --mame` stops after splicing and does not generate expression, PCA,
or plots. A successful standalone `--design-only` run produces only
`design_matrix.tsv`, `contrasts.tsv`, `design.json`, and `glmmTMB.log`;
it does not fit events or generate the remaining regression outputs.

### Regression files in `results/splicing/`

All files in this table are written by a successful full regression run, except
the explicitly optional Excel workbook. Some tables can contain headers only
when there are no applicable results.

| File | Contents and use |
| --- | --- |
| `event_statistics.tsv` | Main results across all event types: one row per event/contrast, event annotations, `p_event`, `q_event`, `status`, component effect summaries, and `Diff events`. Start here for downstream analysis; select the intended contrast before comparing events. |
| `PSI_<TYPE>.txt` | The same event/contrast results split by event type. Despite the historical filename, these are regression results, not the legacy two-group Fisher schema. Sample `*_PSI` columns are added when `--individual-psi` is enabled. |
| `component_statistics.tsv` | One row per tested component/contrast: logit `estimate`, Wald `se`/`ci_low`/`ci_high`, LRT `lr` and `p_component`, dispersion `phi`/`rho`, full/null log likelihoods, eligible sample count, convergence/Hessian flags, `status`, and `diagnostic`. Coverage-filtered events are not fitted. |
| `components.tsv` | Definitions for all input events, including filtered events: component IDs and inclusion/exclusion junction IDs. Use to link a component result back to its junction comparison. |
| `event_filter.tsv` | One row per event/filter group: `n_samples`, `n_eligible`, `n_required`, `group_pass`, and overall `filter_pass`. `filter_group` identifies the categorical cell. Use this table to explain coverage exclusions. |
| `component_predictions.tsv` | Component-level predictions: `kind` distinguishes profiles from effects, `name` identifies the profile/effect, and `value` contains fitted PSI or a PSI-scale contrast. Header-only when no predictions are available. |
| `model_predictions.tsv` | Event-level equal-weight component means for each `kind`/`name`, with `component_min`, `component_max`, and `status`. Incomplete component sets have missing numeric values and `incomplete_components` status; the file itself is still created. Ranges are not confidence intervals. |
| `observed_PSI_<TYPE>.tsv` | Historical observed PSI and junction-count columns, separately from model predictions. Includes input events regardless of the regression coverage filter. |
| `PSI_matrix_sample.txt` | Observed event-by-sample PSI matrix used for PCA; not fitted PSI. This mode does not generate `PSI_matrix_group.txt`. |
| `design_matrix.tsv` | Actual numeric R model matrix, with sample IDs and one column per coefficient. |
| `contrasts.tsv` | One row per coefficient and one column per named contrast, giving its weights. |
| `design.json` | Formula, coefficient names, factor levels, prediction definitions, target hypothesis/population, and R/glmmTMB versions. |
| `sample_metadata.tsv` | Metadata used by the analysis, aligned to the analysis samples. Factor coding and reference levels should be checked in `design.json` and `design_matrix.tsv`. |
| `analysis.json` | Analysis manifest: method/backend, formula, contrasts, FDR method/family/threshold, coverage settings and groups, event/sample counts, effect-filter settings, process count, batch limit, and cache size. Also identifies the regression report mode. |
| `glmmTMB.log` | R backend output from design construction and fitting; inspect alongside component diagnostics if a run or fit fails. |
| `summary.txt` | Counts (`Number`) grouped by `event_type`, `contrast`, `status`, and `Diff events`. Counts across contrasts can include the same event more than once. |
| `factorial_results.xlsx` | Optional (`--excel` / `excel: True`): eight event-type sheets plus `model_predictions`. It does not replace the component, filter, or diagnostic TSV files. |
| `stage_timings.tsv` | Python orchestration stages and elapsed wall-clock seconds for each stage. |
| `component_timings.tsv` | Per-component elapsed/full/null fitting time, optimizer attempts, retries, worker PID, and cache-hit information. |
| `batch_timings.tsv` | Per-batch component/fitted-component counts, cache hits, optimizer attempts, retries, and elapsed seconds. |
| `backend_performance.json` | Backend process/scheduling/cache settings and aggregate workload, cache, optimizer, retry, and elapsed-time measurements. |

Timing files describe different levels of the same computation; do not add
worker times to batch or stage times to estimate total wall time. When no events
pass coverage, fitting/prediction tables and per-component/per-batch timings are
empty, while filter results and backend totals still document the outcome.
Internal `.factorial-*/` batch TSVs, configuration JSON, and `design.rds` are
temporary working files, not retained analysis products.

### Figures and HTML in `plots/`

| File or pattern | Contents and generation condition |
| --- | --- |
| `summary.html` | Self-contained interactive report: command provenance, analysis settings, discovery counts, contrast/status summaries, available PSI effects/PCA, and previews of the first 200 event/prediction rows. Full tables remain in the TSV files. |
| `figure_style.json` | Resolved static font paths, export DPI, and export formats. |
| `png/<name>.png` | Static 300-dpi figure exports. |
| `pdf/<name>.pdf`, `svg/<name>.svg` | Static publication exports; dense scatter points can be rasterized. |
| `interactive/<name>.html` | Plot fragments assembled into the main HTML. These are not standalone pages and do not need to accompany `summary.html`. |
| `data/pca_TPM.html`, `data/pca_PSI.html` | Standalone interactive PCA pages, created only when the corresponding PCA data and experiment table are supplied. |

The figure `<name>` values are:

- `barplot_splicing_summary`: significant-event counts per contrast, stacked by event type.
- `contrast_1`, `contrast_2`, …: conservative logit component effect versus
  −log10(event Q), one figure for each contrast present in the event table.
  Numbering follows contrast appearance in that table; the HTML heading/caption
  gives the contrast name.
- `psi_effect_1`, `psi_effect_2`, …: at most 20 events with the largest absolute
  finite model PSI effect for each available effect name. Numbering follows
  effect-name appearance; figures are omitted when there are no finite effects.
- `pca_TPM` and `pca_PSI`: optional sample PCA figures, requiring their PCA
  tables and an experiment table.

### Other full-pipeline outputs

These files are shared with the bulk workflow; their availability depends on
the selected steps and configuration.

| Path relative to `workdir` | Contents |
| --- | --- |
| `annotation/assembled_annotation.gtf` | Assembled annotation when unannotated-event discovery is enabled, with assembly intermediates in the same directory. |
| `events/EVENT_<TYPE>.txt` | Event definitions used for counting and splicing analysis. |
| `junctions/junctions.bed` | Combined sample count table used by the regression. Counting intermediates vary with workflow and cleanup settings. |
| `results/expression/counts.txt` | Gene count matrix. |
| `results/expression/TPM.txt`, `results/expression/CPM.txt` | Gene abundance matrices; TPM supplies expression PCA. |
| `results/expression/DEG.txt` | Pairwise DESeq2 results only when both expression comparison groups are supplied. This is not factorial splicing output. |
| `results/expression/TPM_CPM.xlsx` | Optional expression workbook when Excel export is enabled. |
| `results/pca/tpm_pca.tsv`, `results/pca/psi_pca.tsv` | Sample principal-component scores for TPM and observed PSI. |
| `results/pca/tpm_contribution.tsv`, `results/pca/psi_contribution.tsv` | Explained variance fractions for each principal component. |
| `report.json` | Pipeline execution provenance, including command, version, input experiment table, and timing. Not generated by standalone `psi.py`. |
| `log/`, `benchmark/` | SnakeShiba rule logs and resource/timing records. These directories are not automatically created by the Python launcher. |

### Reading results and reports

An event has no single fitted regression coefficient: `estimate_min/max` and
`conservative_component_effect` summarize component coefficients, while the
actual estimates remain in the component table. `conservative_component_effect`
is the signed estimate with smallest absolute magnitude, used for the regression
report's horizontal axis. Legacy output columns are unchanged in legacy mode;
consumers of the new mode must use its explicit schema.

The full pipeline's report recognizes `analysis.json` and displays contrasts
without assuming a single reference/alternative PSI pair. It writes the existing
`plots/summary.html` and summary plot paths. The standalone regression report can
also be generated with `src/plots.py` using a parent results directory containing
`splicing/`.

Both Fisher and beta-binomial reports use the same interactive Plotly theme:
white backgrounds, thin axes, restrained colors, and compact typography.
HTML uses larger screen-readable labels and markers than the static publication exports.
The analysis command passed by the pipeline is shown at the top; when regenerating
reports, `--shiba-command` overrides the command recorded in the parent
`report.json`. Missing provenance is explicitly marked as not recorded.
Hover to inspect samples/events, zoom or pan, and click legend entries to
show/hide traces. Plotly JavaScript and figure data are embedded in one
`summary.html`, so no network connection is required. Figures are rendered
as they enter the viewport to bound WebGL resource use; axis/legend state is
retained when scrolling. Static PNG (300 dpi), PDF, and SVG exports remain in
`png/`, `pdf/`, and `svg/`. Browser toolbar export saves the current view as SVG
(WebGL scatter points may be rasterized).

Arial is preferred; when unavailable, static exports use DejaVu Sans and the
browser uses its default sans-serif font. `SHIBA_FONT_DIR` optionally supplies
local font files. `figure_style.json` records font paths for regression exports.

Each contrast has a significance/effect scatter and status counts. Optional
model PSI effects display at most 20 events with the largest absolute effect,
without a significance filter; this selection is stated in the caption.
No additional tests or confidence intervals are inferred. Missing coordinates
are omitted with counts reported, and zero q-values are clipped only for log
display. When PCA files and an experiment table are available, sample PCA
figures are included. Embedded result tables show the first 200 rows; complete
TSVs remain in the analysis output directory.

```bash
python src/plots.py -i /path/to/results -o /path/to/plots \
  -e /path/to/Shiba_experiment.tsv
```

The HTML can be copied independently of its export directories. Large outputs
and validation previews should remain outside version control.

## Shiba / SnakeShiba configuration

Run the complete bulk pipeline with one command:

```bash
python shiba.py config.yaml -p 8
```

The pipeline runs transcript assembly, event extraction, junction counting,
factorial splicing analysis, expression abundance, PCA, and reporting. Use
`--mame` to stop after splicing, or `--start-step 4` to reuse existing events
and junction counts under the configured workdir. R/glmmTMB and jsonlite must
be available to the configured `rscript`; see the installation guide.

Set `stat_method: beta_binomial` to enable factorial regression. Set
`stat_method: fisher` (the default) to select Fisher testing instead. This
single setting controls which splicing analysis is run; no additional Boolean
switch is needed. Supply `workdir`, `gtf`, and `experiment_table` as usual.
The experiment table must contain `sample`, `bam`, and `group`; sample IDs must
match the metadata and count columns. A separate `sample_metadata` file is
optional if the experiment table already contains all formula variables.
Paths are interpreted relative to the invocation directory; absolute paths
are recommended for both Python and Snakemake workflows.

Add these settings to the usual bulk configuration (see the example file):

```yaml
stat_method: beta_binomial
sample_metadata: /path/to/samples.tsv
formula: genotype * treatment + batch
coef:
  - genotypeKO:treatmentDrug
reference_levels:
  genotype: WT
  treatment: Control
categorical: [batch]
minimum_reads: 10
min_sample_fraction: 0.5
# Optional override; default includes every categorical factor, including batch:
# filter_group: [genotype, treatment]
p_adjust: BH
```

`reference_group` and `alternative_group` are not required for regression.
The pipeline disables `ttest` in regression mode; explicitly combining `-t`
with regression at the `psi.py` CLI is rejected. Delete the removed
`beta_binomial` key from existing configurations, even when its value is false. To additionally
run the existing pairwise **gene-expression** analysis, set both
`expression_reference_group` and `expression_alternative_group`. Otherwise
expression abundance and PCA are produced without DESeq2 testing. The splicing
formula is not applied to expression analysis.
