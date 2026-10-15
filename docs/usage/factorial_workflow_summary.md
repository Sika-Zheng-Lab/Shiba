# Title:

A Unified Configuration-Driven Workflow for Fisher and Factorial Beta-Binomial Analysis in Shiba

# Purpose:

This implementation simplifies Shiba's statistical interface by removing the supplementary legacy Python beta-binomial test and making the event-centric R/glmmTMB framework the supported bulk beta-binomial workflow. It also makes the workflow easier to configure and reproduce through consistent method names, matching input templates, and a dedicated MkDocs quick start.

The current work concerns workflow integration, interface consistency, and validation. It does not represent a new biological analysis of the production dataset or a new runtime/memory benchmark.

# Methods:

## Statistical method selection and migration

The configuration and standalone PSI CLI now accept two method values:

- `fisher`: the default when `stat_method` is omitted; retains pairwise Fisher testing and the optional Welch t-test.
- `beta_binomial`: event-centric factorial regression using R/glmmTMB, with biological replicates, fixed-effect formulas, covariates, and scalar linear contrasts.

The old `beta_binomial: true/false` configuration key, `-b/--beta-binomial` PSI options, Python supplementary model, and its dedicated tests were removed. This removal also applies to the single-cell launchers, which shared the old implementation. Factorial regression remains a bulk biological-replicate analysis; this change does not add factorial inference to single-cell mode.

Obsolete configuration keys now cause an explicit migration error. Previous method values `legacy` and `beta-binomial` must be changed to `fisher` and `beta_binomial`. Existing result manifests using the older beta-binomial spelling remain readable by the reporting entry point.

## Pipeline integration

`shiba.py` uses the configuration to route the splicing stage to the selected method. It applies optional defaults before accessing them and uses its own Python interpreter for subprocesses, reducing environment mismatches. The same factorial option mapping is used by the Python and Snakemake launchers.

The full pipeline connects assembly, event extraction, junction counting, splicing analysis, expression abundance, PCA, and reporting. Splicing-only execution and restarting from existing events/counts remain available. Factorial mode disables the supplementary Welch test. Gene-expression differential testing is a separate pairwise operation controlled by the expression group options; it does not use the splicing formula.

## Example inputs and documentation

The supplied factorial configuration uses `genotype * treatment + batch`. Matching TSV templates contain 16 illustrative samples: two genotypes, two treatments, two batches, and two independent biological replicates per cell. The experiment table and metadata share identical unique sample IDs. BAM paths are placeholders; the templates do not include sequencing data or a reference GTF.

The default coefficient is `genotypeKO:treatmentDrug`, with WT and Control as reference levels. Optional JSON templates define the interaction and the drug effects within WT and KO, together with fitted PSI profiles and PSI-scale effects. JSON input is not required for a coefficient-based analysis: `coef` in YAML is sufficient. General linear combinations use `contrast_file`, and explicit prediction definitions use `prediction_grid`.

Coverage filtering remains at least 10 total reads for all required components in the same half of samples, rounded up, within every categorical design cell. With the example metadata, the default includes batch: eight cells with two samples each require at least one qualifying sample per cell. An explicit `filter_group: [genotype, treatment]` instead gives four groups of four samples and requires at least two qualifying samples per group. This changes event eligibility without changing the regression formula.

Execution instructions were moved out of the README into the MkDocs [pipeline command guide](../quickstart/running.md) and [factorial quick start](../quickstart/beta_binomial.md). The README now provides links, and the MkDocs navigation exposes both pages.

## Validation strategy

Validation covered omitted and explicit method selection, rejection of obsolete keys and CLI flags, propagation of formula and grouping options, and construction of all seven pipeline commands. The all-stage command test uses mocked subprocess execution; it does not rerun BAM processing.

A separate integration test executes `shiba.py` from a YAML configuration, starts at the splicing stage with synthetic junction/event inputs, and runs the real R/glmmTMB backend with two workers. It checks successful interaction detection and the new output manifest method value.

The example TSVs were checked for matching unique IDs and balanced cells. Both the coefficient-based example and the optional contrast/prediction JSONs were validated using the actual R design builder. MkDocs navigation targets, local links in the moved pages, and code fences were checked. A full website build was not performed because MkDocs is not installed in the current environment.

# Results:

The final full suite completed 235 tests in 93.121 seconds: **232 passed and three were skipped** because optional dependencies were unavailable. This includes the real R/glmmTMB launcher integration. The command used was:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 R_LIBS_USER=/tmp/shiba-r-library \
/tmp/shiba-factorial-python/bin/python -m unittest discover -s test
```

These Python and R library paths identify the validation environment; users should substitute their own configured environment. `git diff --check` also passed.

The YAML-to-launcher-to-glmmTMB integration successfully detected the simulated positive interaction and produced `analysis.json` with `stat_method: beta_binomial`. Fisher routing worked both with an explicit setting and with no method setting. Removed options failed explicitly rather than silently selecting an unintended analysis.

Both example design checks succeeded. R produced the expected coefficient names: `(Intercept)`, `genotypeKO`, `treatmentDrug`, `batchB2`, and `genotypeKO:treatmentDrug`. The templates therefore provide a consistent starting point for users replacing the illustrative data with their own samples.

## Interpretation

The interaction coefficient tests whether the drug effect differs between KO and WT on the log-odds scale. `treatmentDrug` represents the drug effect within WT; adding the interaction coefficient gives the drug effect within KO. A significant effect in one genotype and a nonsignificant effect in the other does not, by itself, establish a significant interaction.

The optional `delta_delta` prediction is a difference of differences on the PSI scale. It is useful for effect-size interpretation but is not numerically equivalent to the log-odds interaction coefficient. Predictions average over the observed batch distribution in these example profiles. Observed Shiba PSI and component-averaged model predictions should also be distinguished.

The successful simulated interaction confirms that configuration and execution are connected correctly for the tested case. It does not establish power or FDR calibration for every experimental design. The example metadata are illustrative and support no biological conclusion. No new production differential-event counts, biological discoveries, or performance improvements are claimed by this integration update.

# Next step:

1. **Validate the full workflow with a small reproducible BAM dataset.** Exercise counting, factorial analysis, expression abundance, PCA, and reporting together. Compare Python and Snakemake outputs under the same software versions and configuration.
2. **Expand statistical simulation coverage.** Measure false-discovery rates, power, interaction sensitivity, and failure rates across group imbalance, low coverage, dispersion, batch imbalance, and sample sizes representative of intended datasets.
3. **Assess filtering sensitivity computationally.** Compare default categorical-cell filtering with explicitly chosen biological groups; quantify retained events, design estimability, and changes in statistical calls before selecting a grouping strategy.
4. **Evaluate real-data robustness.** Examine component direction consistency, model convergence, fitted versus observed PSI, sample influence, covariate adjustments, and effect-size thresholds for candidate events.
5. **Automate documentation and template checks.** Build MkDocs in CI and validate example YAML/TSV/JSON files against the current parser and R design builder, reducing drift between documentation and implementation.
6. **Revisit resource scaling when workload changes.** Benchmark worker counts and memory for larger sample sets or additional contrasts; this interface update alone does not provide new resource estimates.
