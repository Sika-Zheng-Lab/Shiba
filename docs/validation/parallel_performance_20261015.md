# Title:

Parallel Fisher Testing and Observed PSI Calculation in Shiba: Statistical Equivalence and Real-Data Memory Validation

# Purpose:

This work extends Shiba's performance improvements beyond factorial beta-binomial regression to pairwise Fisher's exact tests and observed PSI calculation. The objectives were to reduce computational overhead, preserve statistical results, and measure whether memory requirements remain practical on the real 65-sample dataset. The validation covers both isolated stages and a complete beta-binomial splicing analysis starting from existing junction counts.

# Methods:

## Implementation

Pairwise Fisher tests now use the requested `-p` process count for sufficiently large event tables. Events are divided into bounded tasks and processed concurrently. Raw results are combined in deterministic order before Benjamini–Hochberg correction is applied once across the original event-type testing family. Chunk boundaries therefore do not change the hypotheses included in FDR correction. Small workloads retain serial execution to avoid startup overhead.

A process-local cache retains up to 10,000 identical 2×2 count-table results. Reuse preserves SciPy's two-sided exact test and odds ratio; it does not merge distinct event hypotheses. MSE, AFE, and ALE significance assignment was vectorized without changing its thresholds or rules.

Factorial observed PSI now honors `-p`. A persistent Python process pool is reused across event types, with the dense junction-count array in shared memory. Workers receive event slices, event ordering is deterministic, and the pool and shared-memory allocation are released even when computation fails. R fitting and Python observed PSI execute in separate stages. These are process-based parallel implementations, rather than Python threading.

The beta-binomial backend already uses persistent R workers, dynamic component scheduling, post-filter event batches, vectorized aggregation, and a bounded exact-result cache. The full-run memory experiment validates these mechanisms together with the new observed PSI implementation.

## Dataset and statistical design

The real dataset contains 65 samples, 3,742,179 junctions, and 159,204 events across SE, FIVE, THREE, MXE, RI, MSE, AFE, and ALE. Its dense int64 junction array alone occupies approximately 1.81 GiB.

The factorial model was `sex + fad * ptbp2ko`, with sex categorical (reference F) and fad and ptbp2ko encoded as numeric binary predictors. Coverage groups were explicitly specified as the combinations of fad and ptbp2ko. This explicit grouping is relevant because the default grouping derives from categorical design factors. Each event had to satisfy the 10-read threshold for every component in the same half of samples, rounded up, within every coverage group. Filtering retained 64,266 events and excluded 94,938 before fitting, leaving 169,768 component records.

Three scalar contrasts were evaluated: interaction (`fad:ptbp2ko`), KO in WT (`ptbp2ko`), and KO in 5xFAD (`ptbp2ko + fad:ptbp2ko`). Factorial FDR correction includes coverage-passing events across event types separately for each contrast; post-filter fitting failures count as p=1 in that family. This differs from the legacy Fisher event-type-specific family, which was preserved.

## Validation and measurement

Tests checked all eight event types, empty or invalid chunks, independent global BH calculations, exact count-table cache results, deterministic PSI output, dispatch, and shared-memory cleanup. The full suite reported 249 passes and three optional-dependency skips (252 tests). After the final persistent-pool refinement, all seven focused tests passed again; the three observed PSI tests also passed with the spawn process-start method.

Single-run synthetic benchmarks compared the original implementation with the optimized implementation. These measurements illustrate possible speedups; they are not replicated estimates of production scaling.

Real-data memory measurements sampled proportional set size (PSS) across the parent and all descendant processes. PSS apportions shared pages rather than counting them in full for every worker. Measurements used fresh processes, included loading, and recorded no memory-access errors. Sampled maxima can miss short peaks and exclude some kernel and filesystem-cache charges.

The Fisher/observed-PSI stage experiments used 1, 4, and 32 workers within a four-CPU, 80-GiB allocation. Consequently, the 32-worker cases measure memory under oversubscription, not 32-CPU throughput. The Fisher experiment compared 9 reference and 10 KO samples in the WT background, testing 129,984 valid group-PSI events. It used the production reader and Fisher functions, but NumPy group-count aggregation and no retained per-sample PSI outputs; it was not the full legacy CLI pipeline.

The complete beta-binomial splicing run used an actual 32-CPU, 80-GiB allocation, native BLAS/OpenMP threads limited to one, a 10,000-entry component cache, and 200 passing events per batch. It included validation, filtering, R/glmmTMB fitting, aggregation, observed PSI, and TSV output. BAM counting, plotting, and Excel export were outside this experiment. The environment used Python 3.12.12, pandas 3.0.3, and R 4.5.2.

# Results:

## Runtime

| Synthetic workload | Original/serial baseline | Optimized, 4 workers | Approximate speedup |
| --- | ---: | ---: | ---: |
| Fisher SE, 8,000 events | 15.879 s | 9.051 s | 1.75× |
| Fisher AFE, 2,000 events | 9.994 s | 7.371 s | 1.36× |
| Observed PSI, 8,000 events × 64 samples | 4.938 s (current serial) | 3.449 s | 1.43× |

Outputs matched in these comparisons. Optimized serial Fisher timings were 16.099 s for SE and 10.021 s for AFE, indicating that the demonstrated benefit primarily came from parallel execution for these workloads.

The complete real-data beta-binomial job finished successfully in **37 minutes 59 seconds**. Model fitting took approximately 1,918 seconds, whereas observed PSI and matrix writing took approximately 25.8 seconds. Fitting therefore remains the main elapsed-time target. The cache reused 27,027 of 169,768 component records (15.9%). This wall time is a measured completion time, not a controlled speedup relative to an earlier run on another node.

## Memory

| Real-data stage experiment | 1 worker | 4 workers | 32 workers |
| --- | ---: | ---: | ---: |
| Fisher stage, peak total PSS | 8.287 GiB | 8.359 GiB | 8.379 GiB |
| Observed PSI stage, peak total PSS | 10.936 GiB | 12.777 GiB | 27.967 GiB |

These stage experiments shared four allocated CPUs. Fisher memory changed little with worker count, whereas observed PSI showed substantial worker-associated overhead despite sharing the dense count array. Shared counts do not eliminate process-local metadata, temporary objects, or retained result tables.

| Complete beta-binomial run, actual 32 CPUs | Peak total process-tree PSS |
| --- | ---: |
| Input and preparation | 11.03 GiB |
| Design stage | 7.23 GiB |
| Model fitting | 9.62 GiB |
| Post-fit processing, including observed PSI and output | **24.18 GiB** |

R alone, including its master and workers, peaked at 6.84 GiB during fitting. Total fitting memory stabilized around 9–10 GiB through the long run and cache turnover, with no substantial progressive growth observed. This supports practical memory use under the tested workload, but does not establish the absence of all memory leaks.

The largest memory demand occurred after fitting. The successful 80-GiB run and measured 24.18-GiB PSS peak support ample headroom in the existing 96-GiB allocation for this splicing workload. A 32-GiB allocation would provide much less margin. These observations are not guarantees for larger datasets, more complex designs, other process-start methods, or the complete BAM-to-report pipeline.

A separate 32-event real-data subset showed fitting peaks of 0.414, 0.918, and 3.997 GiB with 1, 4, and 32 workers. These small-subset measurements demonstrate worker overhead; they must not be interpreted as full-dataset memory scaling or cache-capacity tests.

## Statistical and biological interpretation

Fisher and observed PSI outputs had identical hashes across worker counts within each real-data stage experiment. Comparisons against pre-change Fisher functions also matched across all event types.

The full beta-binomial run matched the existing production results at relative tolerance 1e-8 and absolute tolerance 1e-10 for:

- Event-level statistics: 477,612 event–contrast rows, comparing component counts, status, p-values, effect summaries, sample counts, q-values, and differential-event labels.
- Model predictions: all columns across 387,716 rows.
- Coverage filtering: all columns across 636,816 event–group rows.

The full event-statistics comparison did not include annotation fields or all observed per-sample PSI columns. Those PSI calculations were checked separately in the stage comparisons and focused tests. The 1-, 4-, and 32-worker subset comparisons also matched.

These findings support preservation of the tested statistical outputs while changing scheduling and result reuse. They do not constitute new biological evidence or an independent validation of model calibration. The three contrasts answer different questions: the two KO contrasts estimate the KO effect within each background, whereas the interaction tests whether those effects differ on the model's log-odds scale. Significance in one background alone does not establish an interaction. Model-based event PSI is an unweighted component summary and should also be distinguished from observed Shiba PSI when interpreting effect magnitude.

# Next step:

1. **Measure full-dataset CPU scaling under matched allocations.** Repeat runs at 1, 4, 8, 16, and 32 allocated CPUs, using repeated measurements and the same environment. Record elapsed time, CPU-hours, PSS, and scheduler/cgroup memory to choose an efficient worker count.
2. **Reduce post-fit memory.** Profile junction-index metadata, process-local objects, retained per-event tables, and output assembly. Evaluate streaming or incremental writing and a separate observed-PSI worker limit, checking numerical and ordering equivalence.
3. **Profile and improve model fitting.** Analyze component complexity, optimizer retries, load imbalance, and cache hit rates. Benchmark batch and cache sizes, especially for AFE/ALE events with many components, while retaining convergence diagnostics and existing inference.
4. **Validate the complete legacy Fisher CLI.** Include its group-count aggregation, sample-level outputs, and final writing. The isolated Fisher measurements do not cover those additional costs.
5. **Extend statistical computational validation.** Use null and alternative simulations reflecting observed sample sizes, coverage, dispersion, and group imbalance to assess FDR, power, interaction sensitivity, and convergence failures. Perform sensitivity analyses for coverage thresholds and effect-size cutoffs.
6. **Interpret candidates computationally.** Inspect observed and fitted PSI, effect direction, uncertainty, and component consistency; assess sample influence and available covariates. Compare within-background effects directly through the interaction contrast, then perform appropriately background-controlled functional enrichment.

## Reproducibility records

Compact measurements and source hashes are stored in [parallel_performance_20261015.json](parallel_performance_20261015.json). The optional Linux monitor [monitor_process_memory.py](../../test/benchmarks/monitor_process_memory.py) requires `psutil` and records process-tree PSS, USS, and RSS. For example, from an appropriately allocated compute job:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python test/benchmarks/monitor_process_memory.py --directory /path/to/new/measurement \
  -- python src/psi.py /path/to/junctions.bed /path/to/events /path/to/new/results \
  --stat-method beta-binomial --sample-metadata /path/to/metadata.tsv \
  --formula 'sex + fad * ptbp2ko' --categorical sex \
  --continuous fad --continuous ptbp2ko --reference-level sex=F \
  --filter-group fad --filter-group ptbp2ko \
  --contrast-file /path/to/contrasts.json --prediction-grid /path/to/predictions.json \
  --individual-psi --component-cache-size 10000 --event-batch-size 200 -p 32
```

Use the configured Rscript executable and R library environment described in the installation guide. The original dataset, metadata, large result tables, and raw measurement logs remain outside Git; the compact records do not make the private dataset reproducible without access to those inputs.
