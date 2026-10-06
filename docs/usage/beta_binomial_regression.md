# Event-centric beta-binomial regression

Shiba and SnakeShiba support an optional `stat_method: beta-binomial` mode for
factorial differential splicing. Python extracts integer junction comparisons;
R/glmmTMB fits each component across biological replicates. The default `legacy`
mode, including the older `beta_binomial: True` two-group supplementary test,
retains its existing behavior. These are separate modes.

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
  --stat-method beta-binomial --sample-metadata samples.tsv \
  --formula 'genotype * treatment + batch' \
  --reference-level genotype=WT --reference-level treatment=Control \
  --design-only
```

Then test the interaction:

```bash
python src/psi.py junctions.bed events results \
  --stat-method beta-binomial --sample-metadata samples.tsv \
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

Add these settings to the usual bulk configuration (see the example file):

```yaml
stat_method: beta-binomial
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
Legacy `ttest` and `beta_binomial` switches are ignored by the pipeline's new
mode; explicitly combining them at the `psi.py` CLI is rejected. To additionally
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

R is started for design validation and once for all bounded count batches.
`-p` parallelizes events on Unix; Windows uses sequential fitting. Each fit uses
one internal thread to avoid nested parallelism. Intermediate counts are removed
after the run; fitting logs remain on failure. Completed result files are moved
into the output directory only after the analysis succeeds.

---

## Legacy supplementary two-group beta-binomial implementation

The following describes only `stat_method: legacy` with `beta_binomial: True`.

### Legacy model details

!!! Under development

	This document is currently under development and may contain inaccuracies or incomplete information. Please refer to the latest version on GitHub for updates.

This document summarizes the statistical theory and the current implementation of beta-binomial regression used in Shiba.

The implementation targets event-wise differential splicing testing between two groups, using read counts directly.

## Why beta-binomial

PSI is a ratio derived from count data, and per-sample uncertainty depends strongly on read depth.

If we model PSI alone, depth information is partially lost. Beta-binomial modeling keeps both:

- success counts $k_i$ (inclusion-side reads)
- total counts $n_i$ (effective total reads)

for sample $i$.

This allows overdispersion beyond a simple binomial model while preserving count-scale uncertainty.

## Data structure for one event

For one splicing event and one sample $i$:

- $x_i \in \{0,1\}$: group indicator (reference $=0$, alternative $=1$)
- $k_i$: success reads
- $n_i$: total reads, with $0 \le k_i \le n_i$

The model is fitted per event using all valid samples from both groups.

### Event-wise definition of success and total counts

Shiba computes event-specific $k_i$ and $n_i$ from junction counts as follows.

Let $a,b,c$ denote event-specific junction counts in each sample.

1. FIVE / THREE:

$$
k = a, \quad n = a + b
$$

2. AFE / ALE:

$$
k = \sum a_j, \quad n = \sum a_j + \sum b_j
$$

3. MXE:

$$
k = a_1 + a_2, \quad n = a_1 + a_2 + b_1 + b_2
$$

4. MSE:

$$
k = \sum a_j, \quad n = \sum a_j + c
$$

where $a_j$ are inclusion junctions and $c$ is exclusion junction count.

5. SE:

$$
k = a + b, \quad n = a + b + 2c
$$

This matches the PSI structure

$$
\psi = \frac{(a+b)/2}{(a+b)/2 + c} = \frac{a+b}{a+b+2c} = \frac{k}{n}.
$$

6. RI:

$$
k = s + e, \quad n = s + e + 2a
$$

where $s,e$ are splice-junction counts at intron start/end, and $a$ is intron-retention junction count.

This also preserves

$$
\psi = \frac{(s+e)/2}{(s+e)/2 + a} = \frac{s+e}{s+e+2a} = \frac{k}{n}.
$$

## Statistical model

For sample $i$:

$$
K_i \mid p_i \sim \mathrm{Binomial}(n_i, p_i), \quad p_i \sim \mathrm{Beta}(\alpha_i, \beta_i).
$$

Marginally:

$$
K_i \sim \mathrm{BetaBinomial}(n_i, \alpha_i, \beta_i).
$$

The PMF is

$$
\Pr(K_i = k_i)
= \binom{n_i}{k_i}
\frac{B(k_i + \alpha_i,\; n_i-k_i+\beta_i)}{B(\alpha_i,\beta_i)}.
$$

### Mean and overdispersion parameterization

Shiba uses a mean-dispersion parameterization:

$$
\mu_i = \mathbb{E}[p_i], \quad \rho \in (0,1).
$$

Regression part:

$$
\mathrm{logit}(\mu_i) = \beta_0 + \beta_1 x_i.
$$

Dispersion reparameterization:

$$
\rho = \sigma(\theta) = \frac{1}{1+e^{-\theta}},
$$

$$
\phi = \frac{1-\rho}{\rho},
$$

$$
\alpha_i = \mu_i\phi, \quad \beta_i = (1-\mu_i)\phi.
$$

When $\rho \to 0$, the model approaches binomial-like behavior (low extra-binomial variation).

### What is $\rho$? (core intuition)

In this model, $\rho$ controls **extra-binomial variability**.

If we had a plain binomial model, the success probability is fixed within sample $i$:

$$
K_i \sim \mathrm{Binomial}(n_i, \mu_i), \quad
\mathrm{Var}(K_i) = n_i\mu_i(1-\mu_i).
$$

In beta-binomial, the probability itself fluctuates across biological/technical contexts:

$$
p_i \sim \mathrm{Beta}(\alpha_i,\beta_i), \quad K_i\mid p_i \sim \mathrm{Binomial}(n_i,p_i).
$$

After marginalization, variance is inflated to

$$
\mathrm{Var}(K_i)
= n_i\mu_i(1-\mu_i)\{1 + (n_i-1)\rho\}.
$$

So $\rho$ is the inflation driver:

1. $\rho=0$: no inflation, exactly binomial variance.
2. $\rho>0$: overdispersion; counts are more spread than binomial.
3. Larger $n_i$ amplifies the impact via $(n_i-1)\rho$.

### Equivalent interpretation as intra-class correlation

$\rho$ can also be interpreted as the correlation between two Bernoulli trials within the same sample/event under the beta-binomial hierarchy.

Roughly speaking:

- small $\rho$: reads are close to conditionally independent given $\mu_i$
- large $\rho$: reads co-move more strongly (latent heterogeneity is strong)

This is why beta-binomial is often preferable for RNA splicing counts: technical and biological heterogeneity creates positive correlation and overdispersion that binomial cannot capture.

### Relationship between $\rho$ and $\phi$

Shiba uses

$$
\phi = \frac{1-\rho}{\rho} \quad \Longleftrightarrow \quad \rho = \frac{1}{\phi+1}.
$$

Hence:

1. Large $\phi$ means concentrated Beta distribution around $\mu_i$ and small overdispersion ($\rho \approx 0$).
2. Small $\phi$ means broad Beta distribution and strong overdispersion (large $\rho$).

### Why optimize $\theta$ instead of $\rho$ directly

Implementation uses

$$
\rho = \sigma(\theta) = \frac{1}{1+e^{-\theta}},
$$

so unconstrained $\theta\in\mathbb{R}$ maps smoothly into valid $\rho\in(0,1)$.
This avoids invalid updates during numerical optimization.

## Likelihood used in implementation

For one event, with samples $i=1,\dots,m$, the log-likelihood is

$$
\ell(\Theta)
= \sum_{i=1}^{m}
\left[
\log\binom{n_i}{k_i}
+ \log B(k_i+\alpha_i, n_i-k_i+\beta_i)
- \log B(\alpha_i,\beta_i)
\right].
$$

In code, this is computed via gamma/beta log functions:

$$
\log\binom{n}{k}
= \log\Gamma(n+1)-\log\Gamma(k+1)-\log\Gamma(n-k+1).
$$

Two models are fitted by numerical optimization:

1. Null model $H_0$: $\beta_1 = 0$
2. Full model $H_1$: $\beta_1$ free

Negative log-likelihoods are minimized with L-BFGS-B.

## Hypothesis test

Primary hypothesis:

$$
H_0: \beta_1 = 0
\quad\text{vs}\quad
H_1: \beta_1 \ne 0.
$$

Likelihood-ratio statistic:

$$
\Lambda = 2\{\ell(\widehat\Theta_1) - \ell(\widehat\Theta_0)\}.
$$

Implementation uses equivalent form with minimized negative log-likelihoods:

$$
\Lambda = 2\{\mathrm{NLL}_{0} - \mathrm{NLL}_{1}\}, \quad \Lambda \leftarrow \max(\Lambda,0).
$$

P-value per event:

$$
p_{\beta} = \Pr\left(\chi^2_{(1)} \ge \Lambda\right).
$$

Then BH correction is applied across valid events to create $q_{\beta}$.

## Practical safeguards in Shiba implementation

The implementation includes several robustification steps.

1. Sample count requirement:

- if either group has fewer than 2 valid samples, result is `NaN`.

2. Early exits for near-constant events:

- if variance of $k_i/n_i$ is extremely small, return $\Lambda=0$.
- if group means are almost identical, return $\Lambda=0$.

3. Parameter clipping and bounds:

- $\mu$ clipped to $(10^{-10}, 1-10^{-10})$
- $\rho$ clipped to $(10^{-8}, 1-10^{-8})$
- optimization bounds:
  - null: $\beta_0\in[-20,20],\;\theta\in[-10,10]$
  - full: $\beta_0\in[-20,20],\;\beta_1\in[-20,20],\;\theta\in[-10,10]$

4. Numerical validity checks:

- if optimizer returns non-finite objective values, event result is `NaN`.

5. P-value underflow handling:

- when $p_{\beta}$ underflows to zero, it is replaced with `np.finfo(float).tiny`.

## Relationship to PSI output

Important distinction:

- Statistical inference now uses direct count inputs $(k,n)$.
- PSI is still reported as an interpretable summary metric.

So PSI remains useful for biological interpretation, while hypothesis testing is count-aware and depth-aware.

## End-to-end flow in Shiba

For each event type:

1. Compute group-level PSI and detect candidate differential events.
2. For candidate events, compute sample-level PSI and count columns (`_success`, `_total_reads`).
3. Run beta-binomial LRT to get `p_beta`.
4. Apply BH correction to get `q_beta`.
