# Output files of scShiba/SnakeScShiba

All paths are relative to the configured `workdir`. Single-cell mode uses
Fisher testing on counts aggregated by cell group; factorial beta-binomial
regression is available only in the bulk workflow.

| Path | Contents |
| --- | --- |
| `events/EVENT_<TYPE>.txt` | Annotated event definitions for SE, FIVE, THREE, MXE, RI, MSE, AFE, and ALE. |
| `junctions/junctions.bed` | Tosa junction and RI-boundary counts aggregated by group. |
| `results/PSI_<TYPE>.txt` | Event annotations, observed group PSI, and differential results. In PSI-only mode, contains observed PSI without tests. |
| `results/summary.txt` | Differential-event counts by event type, annotation, and direction; differential runs only. |
| `results/PSI_matrix_sample.txt` | PSI matrix written only with `only_psi: True`; columns represent the groups in the aggregated count table despite the historical filename. |
| `results/results.xlsx` | Optional Excel export. |
| `plots/summary.html` | Self-contained interactive report with command provenance, summary counts, and volcano/PSI scatter plots for all eight types. PSI-only runs instead show labeled table previews. |
| `plots/png/`, `plots/pdf/`, `plots/svg/` | Static figure exports for differential runs. |
| `plots/interactive/` | Plot fragments embedded in the main HTML. |
| `plots/data/` | Standalone event scatter and volcano HTML pages. |
| `report.json` | Pipeline command, version, input experiment table, and execution timing. |
| `benchmark/`, `log/` | Rule resource measurements and logs added by SnakeScShiba. |

Results are saved directly in `results/`, **not** `results/splicing/`.
There are no expression or PCA outputs. Intermediate Tosa files may be removed
by workflow cleanup.

See the [Fisher column reference](shiba.md#column-description-of-psi_txt),
[single-cell quick start](../quickstart/diff_splicing_sc.md), and
[report regeneration instructions](../usage/scshiba.md#interactive-summary-report).
