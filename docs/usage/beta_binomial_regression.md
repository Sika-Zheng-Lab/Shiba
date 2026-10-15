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

| File | Contents |
| --- | --- |
| `PSI_<TYPE>.txt`, `event_statistics.tsv` | Annotations and one row per event/contrast; P/Q, status, component effect range and calls |
| `component_statistics.tsv` | Contrast estimates, Wald SE/CI, LRT, dispersion, sample count and fitting diagnostics |
| `components.tsv` | Junction IDs and component definitions |
| `event_filter.tsv` | All events: group definitions, eligible/total/required sample counts, per-group pass and overall filter pass |
| `component_predictions.tsv` | Each component's profile PSI and PSI effects |
| `model_predictions.tsv` | Event mean predictions/effects and component ranges; missing if incomplete |
| `observed_PSI_<TYPE>.tsv` | Historical observed PSI and junction counts |
| `PSI_matrix_sample.txt` | Observed PSI matrix for PCA |
| `design_matrix.tsv`, `contrasts.tsv`, `design.json` | Actual design, contrast weights, factor levels, prediction definitions and versions |
| `sample_metadata.tsv`, `analysis.json`, `glmmTMB.log` | Reproducibility metadata, analysis settings and backend log |
| `summary.txt` | Counts by event type, contrast, status and significance |
| `factorial_results.xlsx` | Optional event and model prediction sheets (`--excel`) |
| `stage_timings.tsv` | Python orchestration stages and elapsed wall-clock seconds |
| `component_timings.tsv` | Per-component fitting time, full/null optimizer attempts, retries, worker PID and cache hits |
| `batch_timings.tsv`, `backend_performance.json` | Batch wall times and backend workload/cache/retry totals |

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

## Validation and limitations

Run the Python tests and real-backend integration tests:

```bash
python -m unittest discover -s test -p 'test_factorial.py' -v
```

Real-backend tests require glmmTMB/jsonlite and are mandatory in the dedicated CI
job. They cover a known interaction, null, contradictory junctions, unusable
counts, general/reversed contrasts, ordering, parallel execution and designs.
`test/simulate_factorial.py` provides a reproducible calibration/benchmark run.

LRT P values are asymptotic. Small biological replicate counts, separation,
near-binomial dispersion or poor coverage can make inference unreliable even
when an optimizer returns success. The backend checks convergence, Hessian and
likelihood ordering and retries failed fits with a second optimizer. It does
not silently switch to a different distribution. This implementation does not
provide dispersion shrinkage or a parametric-bootstrap P value. Calibration on
the intended sample sizes and depths is essential before scientific use.

A preliminary run on 2026-10-05 with R 4.5.2 / glmmTMB 1.1.13 used seed 731,
200 SE events (100 interaction-null and 100 non-null), depth 100 and true rho
0.05. The two inclusion components shared identical counts. Results were:

| Replicates per cell | Null P < .05 | Realized false discovery proportion after BH | Power | Median estimated rho |
| --- | --- | --- | --- | --- |
| 8 | 9/100 (95% interval 4.2–16.4%) | 3.0% | 97% | 0.0406 |
| 16 | 5/100 (95% interval 1.6–11.3%) | 2.0% | 100% | 0.0451 |

These are two individual simulated datasets, not repeated-simulation estimates
of expected FDR or guarantees for real RNA-seq experiments. The smaller-sample
run illustrates why asymptotic calibration and dispersion bias need attention.
Different coverage, junction dependence and replicate counts require further
calibration. Run both examples with `test/simulate_factorial.py --events 200
--replicates 8` (or `16`), `--seed 731`, and a new `--output` directory.

### Parallel execution and performance

Use `-p N` with `shiba.py` or `src/psi.py`, or `--cores N` with SnakeShiba.
The default for the Python commands is one process. R is started for design
validation and once for fitting. The fitting stage keeps workers alive across
batches and dynamically assigns component fits to them, so a large AFE/ALE
event does not occupy just one worker. Unix uses a fork cluster; Windows uses
a PSOCK cluster. Each component's full model and null models run sequentially
within one worker, with one internal fitting thread.

Observed PSI calculation also uses `-p N`, reusing one Python worker pool
across event types and sending each worker only its event slice. Junction counts are shared between
Python workers through shared memory, which is released when the calculation
finishes or fails. Event output order is preserved; empty event types do not
start a worker pool. R fitting and observed PSI run in separate stages, so their
worker pools do not run concurrently.

The event-wide eligible sample mask is established **before** components are
scheduled. Parallel scheduling does not change that mask, the likelihood,
component consistency checks, or the FDR family. Passing events are packed
across input batches and event types into fitting batches of up to
`max(event_batch_size, 4 * processes)` events. The default `event_batch_size`
is 200; it also bounds the number of input events extracted at a time. Memory
depends on component counts and sample counts as well as the number of events.

Within a run, exact component results are reused when the ordered sample mask
and included inclusion/exclusion counts match. The design, contrasts,
prediction population, and optimizer settings are fixed for the run. Keys use
exact serialization, without rounding or approximate matching. Results are
relabeled for each event/component, which remains a separate entry in testing
and FDR correction. Diagnostics, including failures, are retained on reuse.
The bounded FIFO cache stores compact results, not fitted model objects.
Use `--component-cache-size 0` to disable both within-batch deduplication and
cross-batch caching for comparisons; the default is 10,000 entries. In YAML:

```yaml
event_batch_size: 200
component_cache_size: 10000
```

`stage_timings.tsv` separates input validation, design, filtering/batching,
backend fitting, reading, event aggregation, prediction aggregation/output,
observed PSI, and final output. `component_timings.tsv` counts actual optimizer
attempts (initial fits plus retries). Cache-hit rows have zero time/attempts and
no worker PID; the original diagnostics remain in `component_statistics.tsv`.
Times summed over workers are not the wall-clock duration of a parallel run.
Use batch wall times and `backend_performance.json` for overall backend timing.
An all-filtered run writes empty timing tables and zero backend workload.

Event and prediction aggregation use vectorized grouped operations, and
component output is written in batches. These changes preserve statistical
definitions; floating-point summaries can differ in their last digits.
Intermediate counts are removed after the run; fitting logs remain on failure.
Completed result files are moved into the output directory only after the
analysis succeeds. Timings describe completed runs; a failed run retains its
backend log rather than publishing partial result tables.

See the [parallel performance and real-data memory validation report](../validation/parallel_performance_20261015.md) for benchmark conditions, statistical equivalence checks, memory measurements, and computational next steps.

A complete [beta-binomial example configuration](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/Shiba_SnakeShiba_beta_binomial_config_example.yaml) is available. `stat_method` accepts `fisher` or `beta_binomial`; omitting it runs Fisher testing.
