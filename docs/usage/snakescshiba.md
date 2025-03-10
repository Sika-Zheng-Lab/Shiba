# SnakeScShiba Usage

``` bash
snakemake -s snakescshiba.smk \
--configfile config.yaml \
--cores 32 \
--use-singularity \
--singularity-args "--bind $HOME:$HOME" \
--rerun-incomplete
```

Please check the [Quick Start](../quickstart/diff_splicing_sc.md/#1-prepare-inputs_1){ data-preview } to learn how to prepare the `config.yaml`.

<figure markdown="span">
	![SnakeScShiba rulegraph](https://github.com/Sika-Zheng-Lab/Shiba/blob/mkdocs/img/SnakeScShiba_rulegraph.svg?raw=true){ width="500" align="center" }
	<figcaption>SnakeScShiba rulegraph</figcaption>
</figure>

This rulegraph was made by [snakevision](https://github.com/OpenOmics/snakevision).
