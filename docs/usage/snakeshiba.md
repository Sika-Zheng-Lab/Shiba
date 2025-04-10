# SnakeShiba Usage

This **snakeshiba.smk** script is written for use with Snakemake, a workflow management system that allows for the creation of reproducible and scalable data analyses. Snakemake uses a Python-based domain-specific language to define rules that specify how to generate output files from input files. These rules are executed in a directed acyclic graph (DAG) to ensure efficient and correct execution.

Key Features of Snakemake:

- Automatically determines the order of execution based on dependencies.
- Supports parallel execution on local machines, clusters, and cloud environments.
- Provides built-in support for logging, benchmarking, and resource management.

!!! Tip

	For more information about Snakemake, visit:

	- Official Documentation: [https://snakemake.readthedocs.io/](https://snakemake.readthedocs.io/)
	- GitHub Repository: [https://github.com/snakemake/snakemake](https://github.com/snakemake/snakemake)
	- Tutorials and Examples: [https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)

You can run the script using the following command:

``` bash
snakemake -s snakeshiba.smk \
--configfile config.yaml \
--cores 4 \
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
