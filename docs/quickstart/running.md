# Running Shiba pipelines

Run these commands from the Shiba checkout after [installing dependencies](../installation.md) and preparing `config.yaml`. See the [bulk](diff_splicing_bulk.md), [single-cell](diff_splicing_sc.md), or [factorial beta-binomial](beta_binomial.md) quick start for input preparation.

Fisher analysis also uses `-p N` to distribute event tests across processes.
BH correction is applied after joining results, across the same event-type
family as before. Small workloads (fewer than 256 input events per type) run
serially; larger workloads use up to `N` workers. Each process caches up to
10,000 identical 2×2 tables. This cache is separate from the factorial
`component_cache_size` setting.

## Shiba

```bash
python ./shiba.py -p 4 config.yaml
```

## MameShiba

A lightweight version of Shiba.

```bash
python ./shiba.py --mame -p 4 config.yaml
```

## SnakeShiba

Snakemake-based workflow of Shiba.

```bash
snakemake -s snakeshiba.smk --configfile config.yaml --cores 4 --use-singularity
```

## scShiba

A single-cell RNA-seq version of Shiba.

```bash
python ./scshiba.py -p 4 config.yaml
```

## SnakeScShiba

Snakemake-based workflow of scShiba.

```bash
snakemake -s snakescshiba.smk --configfile config.yaml --cores 4 --use-singularity
```
