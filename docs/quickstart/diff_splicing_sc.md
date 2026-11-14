# Differential RNA splicing analysis with single-cell/nucleus RNA-seq data

---

## Before you start

- Map sc(sn)RNA-seq reads to the reference genome and keep a coordinate-sorted BAM/CRAM with `CB` (cell barcode) and `UB` (UMI) tags. Tosa uses `CB` to select cells and `UB` to deduplicate molecules.
- Download a gene annotation file of your interest in GTF format.

---

## Installation

Follow the [installation guide](../installation.md) and activate the environment
containing Shiba. The commands below use the Conda-installed launchers. For a
source checkout, use the [source invocation](../installation.md#locating-scripts-and-workflows)
instead. These examples require the v1.0.0 implementation, including Tosa.

---

## scShiba

### 1. Prepare inputs

`experiment.tsv`: A **tab-separated** table with one row per library. `alignment` points to a BAM or CRAM file; `barcode` points to the library's barcode/group table.

``` text
sample alignment barcode
run1 /path/to/run1/Aligned.sortedByCoord.out.bam /path/to/barcodes_run1.tsv
run2 /path/to/run2/Aligned.sortedByCoord.out.bam /path/to/barcodes_run2.tsv
```

`barcodes.tsv` is a **tab-separated** text file of barcode and group name like this. Barcodes must match the `CB` tags, including any suffix such as `-1`:

``` text
barcode group
TTTGTTGTCCACACCT Cluster-1
TCAAGACCACTACAGT Cluster-1
TATTTCGGTACAGTAA Cluster-1
ATCCTATGTTAATCGC Cluster-1
ATCGATGAGTTTCTTC Cluster-2
ATCGATGGTCTTGCTC Cluster-2
TATGTTCGTCAGGCAA Cluster-2
ATCGCCTAGACTCGAG Cluster-2
...
```

!!! warning "Make sure to use tabs"

    If you copy and paste the above example, your experiment.tsv file may contain **spaces** instead of tabs, which will causes an error when you run **scShiba**. Please make sure that you are using a **tab** character between the columns.

`config.yaml`: A yaml file of the configuration.

This guide describes the Tosa-enabled v1.0.0 implementation. The v0.8.2 package
uses the earlier STARsolo `SJ` input. To try generated data from a source
checkout, run `python test/make_sc_tosa_fixture.py /path/to/fixture` from the
checkout root and use its `config.yaml`. This test helper is not packaged.

``` yaml
workdir:
  /path/to/workdir # (1)!
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf # (2)!
experiment_table:
  /path/to/experiment.tsv # (3)!

# Tosa counting
minimum_anchor_length: 8
boundary_anchor_length: 6
minimum_intron_length: 20
maximum_intron_length: 500000
strand: unstranded

# PSI calculation
only_psi:
  False # (4)!
fdr:
  0.05 # (5)!
delta_psi:
  0.1 # (6)!
reference_group:
  Cluster-1 # (7)!
alternative_group:
  Cluster-2 # (8)!
minimum_reads:
  10 # (9)!
excel:
  False # (10)!
```

1. The working directory where the output files will be saved. Please make sure that you have write permission to this directory.
2. The path to the gene annotation file in GTF format.
3. The path to the `experiment.tsv` file.
4. Set to `True` if you want to skip the differential analysis and only calculate PSI values for each sample.
5. Significance threshold for differential splicing analysis.
6. Minimum difference in PSI values between groups to be considered significant.
7. Reference group for differential splicing analysis.
8. Alternative group for differential splicing analysis.
9. Minimum number of reads required to calculate PSI values.
10. Set to `True` if you want to generate a file of splicing analysis results in excel format.

### 2. Run

``` bash
scshiba.py -p 4 config.yaml
```

You are going to use 4 threads for parallelization. You can change the number of threads by changing the `-p` option.

Results now include `PSI_RI.txt` for intron-retention analysis, alongside the existing seven event types.

!!! bug "Did you encounter any problems?"

	You can run **scShiba** with the `--verbose` option to see the debug log. This will help you to find the problem.
	```bash
	scshiba.py --verbose -p 4 config.yaml
	```
	If you continue to encounter issues, please don't hesitate to [open an issue](https://github.com/Sika-Zheng-Lab/Shiba/issues) on GitHub. The community and developers are here to help!

---

## SnakeScShiba

A snakemake-based workflow of **scShiba**. This is useful for running **scShiba** on a cluster. Snakemake automatically parallelizes the jobs and manages the dependencies between them.

### 1. Prepare inputs

`experiment.tsv`: The same `sample`, `alignment`, `barcode` table used by **scShiba**.

`config.yaml`: A yaml file of the configuration. This is the same as the input for **scShiba** but with the addition of the `container` field.

The container reference below must resolve to a published release matching your source. Until that release is available, run Snakemake locally with the source environment or build `docker/Dockerfile_develop` and set `container` to your matching SIF image.

``` yaml
workdir:
  /path/to/workdir # (1)!
container: # This field is required for SnakeScShiba
  docker://naotokubota/shiba:v1.0.0 # (2)!
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf # (3)!
experiment_table:
  /path/to/experiment.tsv # (4)!

# Tosa counting
minimum_anchor_length: 8
boundary_anchor_length: 6
minimum_intron_length: 20
maximum_intron_length: 500000
strand: unstranded

# PSI calculation
only_psi:
  False # (5)!
fdr:
  0.05 # (6)!
delta_psi:
  0.1 # (7)!
reference_group:
  Cluster-1 # (8)!
alternative_group:
  Cluster-2 # (9)!
minimum_reads:
  10 # (10)!
excel:
  False # (11)!
```

1. The working directory where the output files will be saved. Please make sure that you have write permission to this directory.
2. The Docker image of **Shiba**.
3. The path to the gene annotation file in GTF format.
4. The path to the `experiment.tsv` file.
5. Set to `True` if you want to skip the differential analysis and only calculate PSI values for each sample.
6. Significance threshold for differential splicing analysis.
7. Minimum difference in PSI values between groups to be considered significant.
8. Reference group for differential splicing analysis.
9. Alternative group for differential splicing analysis.
10. Minimum number of reads required to calculate PSI values.
11. Set to `True` if you want to generate a file of splicing analysis results in excel format.

### 2. Run

Install Snakemake and Apptainer/Singularity, set
[`SHIBA_HOME`](../installation.md#locating-scripts-and-workflows), and select a
matching container image as described in the [installation guide](../installation.md#snakemake).
For execution in the activated Conda environment, omit `--use-singularity`.

``` bash
snakemake -s "$SHIBA_HOME/snakescshiba.smk" \
--configfile config.yaml \
--cores 16 \
--use-singularity \
--rerun-incomplete
```
