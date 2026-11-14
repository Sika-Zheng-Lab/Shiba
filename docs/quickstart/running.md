# Running Shiba pipelines

After [installing Shiba](../installation.md), activate its Conda environment and
prepare `config.yaml`. These examples use the installed launchers. For a source
checkout, substitute the [source invocation](../installation.md#locating-scripts-and-workflows).
Set `SHIBA_HOME` as described there before running a Snakemake workflow. See the [bulk](diff_splicing_bulk.md), [single-cell](diff_splicing_sc.md), or [factorial beta-binomial](beta_binomial.md) quick start for input preparation.

Fisher analysis also uses `-p N` to distribute event tests across processes.
BH correction is applied after joining results, across the same event-type
family as before. Small workloads (fewer than 256 input events per type) run
serially; larger workloads use up to `N` workers. Each process caches up to
10,000 identical 2×2 tables. This cache is separate from the factorial
`component_cache_size` setting.

RI boundary counting defaults to `boundary_anchor_length: 6` (six aligned bases
on each side). Explicit config values override this default. When changing it,
regenerate Tosa counts from the alignments with the same setting for counting
and merging; do not reuse intermediate boundary tables made with another value.

## Shiba

```bash
shiba.py -p 4 config.yaml
```

## MameShiba

A lightweight version of Shiba.

```bash
shiba.py --mame -p 4 config.yaml
```

## SnakeShiba

Snakemake-based workflow of Shiba. The example uses the active Conda environment.
For container execution, follow the [container setup](../installation.md#snakemake)
and add `--use-apptainer`.

```bash
snakemake -s "$SHIBA_HOME/snakeshiba.smk" --configfile config.yaml --cores 4
```

## scShiba

A single-cell RNA-seq version of Shiba.

```bash
scshiba.py -p 4 config.yaml
```

## SnakeScShiba

Snakemake-based workflow of scShiba.

```bash
snakemake -s "$SHIBA_HOME/snakescshiba.smk" --configfile config.yaml --cores 4
```
