# SnakeShiba Usage

``` bash
snakemake -s snakeshiba.smk \
--configfile config.yaml \
--cores 32 \
--use-singularity \
--singularity-args "--bind $HOME:$HOME" \
--rerun-incomplete
```

Please check the [Quick Start](../quickstart/diff_splicing_bulk.md/#1-prepare-inputs_1) to learn how to prepare the `config.yaml`.

<figure markdown="span">
	![SnakeShiba rulegraph](https://github.com/Sika-Zheng-Lab/Shiba/blob/mkdocs/img/SnakeShiba_rulegraph.svg?raw=true){ width="500" align="center" }
	<figcaption>SnakeShiba rulegraph</figcaption>
</figure>

This rulegraph was made by [snakevision](https://github.com/OpenOmics/snakevision).
