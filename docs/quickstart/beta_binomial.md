# Factorial beta-binomial analysis

Use the [dedicated config](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/Shiba_SnakeShiba_beta_binomial_config_example.yaml)
with these matching templates:

| File | Purpose |
| --- | --- |
| [samples.tsv](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/beta_binomial/samples.tsv) | Sample metadata for `genotype * treatment + batch` |
| [Shiba_experiment.tsv](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/beta_binomial/Shiba_experiment.tsv) | Matching sample IDs, BAM paths, and pipeline group labels |
| [contrasts.json](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/beta_binomial/contrasts.json) | Optional interaction and drug effects within WT and KO |
| [predictions.json](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/example/beta_binomial/predictions.json) | Optional fitted PSI profiles and drug/interaction PSI effects |

These are **input templates, not a runnable sequencing dataset**. They contain
16 illustrative biological replicates: WT/KO × Control/Drug × B1/B2, with two
replicates per cell. Replace the rows with your independent biological samples;
do not duplicate samples to fill the design. BAM files and a matching reference
GTF must be supplied separately. Install the dependencies in the
[installation guide](../installation.md), including R/glmmTMB and jsonlite.

From the Shiba checkout, prepare a directory (choose a new directory to avoid
overwriting existing configuration files):

```bash
mkdir -p factorial_example
cp example/Shiba_SnakeShiba_beta_binomial_config_example.yaml factorial_example/config.yaml
cp example/beta_binomial/samples.tsv factorial_example/samples.tsv
cp example/beta_binomial/Shiba_experiment.tsv factorial_example/Shiba_experiment.tsv
cp example/beta_binomial/contrasts.json factorial_example/contrasts.json
cp example/beta_binomial/predictions.json factorial_example/predictions.json
```

Before running, edit the copied files:

1. Set absolute `workdir` (output directory), `gtf`, `experiment_table`, and
   `sample_metadata` paths in `config.yaml`, and actual BAM paths in
   `Shiba_experiment.tsv`. Use the same unique sample IDs in both TSV files.
2. Adjust `formula`, `categorical`, `reference_levels`, and `coef` to your
   experimental design. The supplied config tests `genotypeKO:treatmentDrug`.
   `group` labels are used by the pipeline; regression uses the metadata factors.
3. Keep `stat_method: beta_binomial`. `stat_method: fisher`, or omitting it,
   selects Fisher testing. The removed top-level `beta_binomial: true/false`
   switch must not be used.

Run the complete pipeline, including expression abundance, PCA, and reports:

```bash
python shiba.py factorial_example/config.yaml -p 8
```

For splicing only, or to reuse existing events and junction counts under the
configured `workdir`, respectively:

```bash
python shiba.py factorial_example/config.yaml -p 8 --mame
python shiba.py factorial_example/config.yaml -p 8 --start-step 4 --mame
```

To run all three example contrasts, comment out `coef` and set `contrast_file`
to the absolute path of the copied `contrasts.json`. Set `prediction_grid` to
the copied `predictions.json` for the corresponding PSI profiles/effects. The
interaction tests a difference of effects on the log-odds scale; `delta_delta`
is a PSI-scale difference of differences and is not the same numerical effect.
Batch is averaged over the observed samples in these prediction profiles.

Coverage filtering defaults to at least 10 reads in the same half of samples
in every categorical design cell, including batch. In the example, each of the
eight cells has two samples, so at least one sample per cell must meet coverage
for all components. To explicitly group coverage by genotype × treatment only,
uncomment `filter_group: [genotype, treatment]`; each of those four groups then
has four samples, requiring at least two. This changes the filter, not the model.

Results are written under `<workdir>/results/splicing/`, including
`event_statistics.tsv`, `event_filter.tsv`, `design.json`, and, when requested,
`model_predictions.tsv`. See the [regression guide](../usage/beta_binomial_regression.md)
for contrast definitions and interpretation. Gene-expression differential testing
is a separate pairwise analysis enabled by `expression_reference_group` and
`expression_alternative_group`; the splicing formula is not applied to expression.
