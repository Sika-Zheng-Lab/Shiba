# SnakeScShiba Usage

This **snakescshiba.smk** script is written for use with Snakemake, a workflow management system that allows for the creation of reproducible and scalable data analyses. Snakemake uses a Python-based domain-specific language to define rules that specify how to generate output files from input files. These rules are executed in a directed acyclic graph (DAG) to ensure efficient and correct execution.

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
snakemake -s snakescshiba.smk \
--configfile config.yaml \
--cores 4 \
--use-singularity \
--singularity-args "--bind $HOME:$HOME" \
--rerun-incomplete
```

Please check the [Quick Start](../quickstart/diff_splicing_sc.md#1-prepare-inputs_1){ data-preview } to learn how to prepare the `config.yaml`.

## Workflow

Event extraction → per-library Tosa counting → group aggregation → PSI/Fisher analysis → interactive report.

Snakemake schedules independent rules according to their dependencies.


## Interactive summary report

The workflow automatically creates `plots/summary.html` in the work directory.
This self-contained HTML uses the same interactive style as Shiba: enlarged text,
muted colors, hover labels, zoom/pan, legend toggling, and Arial with browser-default
font fallback. It includes the analysis command, differential-event counts, and
volcano and PSI scatter plots for all eight event types, including RI.
PNG, PDF, and SVG exports are saved alongside the report.

Single-cell PSI values are aggregated by cell group. This report does not add
biological replication or change the Fisher tests. Bulk TPM/PSI PCA panels are
not included. With `only_psi: True`, the report instead shows up to 200 rows per
event type and explicitly states that no differential testing was performed.

To regenerate the report from existing results:

```bash
python src/plots.py --single-cell -i /path/to/workdir/results -o /path/to/workdir/plots
```

Add `--only-psi` for PSI-only results. The command is retrieved from the saved
`report.json`; use `-s 'original analysis command'` to supply it explicitly.
