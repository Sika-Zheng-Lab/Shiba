# Differential RNA splicing analysis with bulk RNA-seq data

---

## Before you start

- Perform mapping of RNA-seq reads to the reference genome and generate bam files with their index files (`.bai`)  by software such as [STAR](https://github.com/alexdobin/STAR) and [HISAT2](https://daehwankimlab.github.io/hisat2/).
    - You can download test RNA-seq bam files with their index (two replicates for reference and alternative groups) mapped by STAR on the mouse genome from [here](https://zenodo.org/records/14976391).
- Download a gene annotataion file of your interest in GTF format.

---

## Installation

- **Shiba**:

``` bash
# Install Shiba with conda
conda create -n shiba -c conda-forge -c bioconda shiba
# Activate the conda environment
conda activate shiba
# Install styleframe for generating outputs in Excel format (optional)
pip install styleframe==4.1
```

- **MameShiba**, a lightweight version of **Shiba**:

``` bash
# Install MameShiba with conda
conda create -n mameshiba -c conda-forge -c bioconda mameshiba
# Activate the conda environment
conda activate mameshiba
```

---

## Shiba

### 1. Prepare inputs

`experiment.tsv`: A **tab-separated** text file of sample ID, path to bam files, and groups for differential analysis.

``` text
sample  bam group
Ref_1 /path/to/workdir/bam/Ref_1.bam  Ref
Ref_2 /path/to/workdir/bam/Ref_2.bam  Ref
Alt_1 /path/to/workdir/bam/Alt_1.bam  Alt
Alt_2 /path/to/workdir/bam/Alt_2.bam  Alt
```

!!! warning "Make sure to use tabs"

    If you copy and paste the above example, your experiment.tsv file may contain **spaces** instead of tabs, which will causes an error when you run **Shiba**. Please make sure that you are using a **tab** character between the columns.

!!! info "Shiba supports long-read RNA-seq data"

    If you have **long-read RNA-seq data** (i.e., PacBio or ONT), please add the 4th column to the `experiment.tsv` file with the value `long` for long-read data and `short` for short-read data. For example:

    ``` text
    sample bam group technology
    Ref_1 /path/to/workdir/bam/Ref_1.bam Ref short
    Ref_2 /path/to/workdir/bam/Ref_2.bam Ref long
    Alt_1 /path/to/workdir/bam/Alt_1.bam Alt short
    Alt_2 /path/to/workdir/bam/Alt_2.bam Alt long
    ```

    The 4th column is optional. If you do not have long-read data, you can omit the 4th column. Blank values are also accepted and will be treated as `short`.

    !!! tip "Use long-read data only for discovery of alternative RNA splicing events"

        If you want to use long-read RNA-seq data only for discovery of alternative RNA splicing events and **NOT** for differential analysis, you can set the 3rd column to different values than that of short-read data. For example, if you want to perform differential splicing analysis between `Ref` and `Alt` groups using short-read data, you can set `Ref` and `Alt` for short-read data, and set `Ref_long` and `Alt_long` for long-read data, so that the long-read data will be used only for transcript assembly and not for differential analysis.

`config.yaml`: A yaml file of the configuration.

``` yaml
workdir:
  /path/to/workdir # (1)!
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf # (2)!
experiment_table:
  /path/to/experiment.tsv # (3)!
unannotated:
  True # (4)!

# Junction read filtering
minimum_anchor_length:
  6 # (5)!
minimum_intron_length:
  70 # (6)!
maximum_intron_length:
  500000 # (7)!
strand:
  XS # (8)!

# PSI calculation
only_psi:
  False # (9)!
only_psi_group:
  False # (10)!
fdr:
  0.05 # (11)!
delta_psi:
  0.1 # (12)!
reference_group:
  Ref # (13)!
alternative_group:
  Alt # (14)!
minimum_reads:
  10 # (15)!
individual_psi:
  True # (16)!
ttest:
  True # (17)!
excel:
  False # (18)!
```

1. The path to the working directory. This is where the output files will be saved. Please make sure that you have write permission to this directory.
2. The path to the gene annotation file in GTF format.
3. The path to the `experiment.tsv` file.
4. True if you want to include unannotated splicing events in the analysis. If False, only annotated events are considered.
5. Junctions having a minimum overlap of this much on both ends are reported.
6. Minimum length of the intron sequence.
7. Maximum length of the intron sequence.
8. Strand specificity of RNA library preparation, where the options XS, use XS tags provided by aligner; RF, first-strand; FR, second-strand.
9. True if you want to skip the differential analysis and only calculate PSI values for each sample.
10. True if you want to skip the differential analysis and only calculate PSI values for each group.
11. Significance threshold for differential splicing analysis.
12. Minimum difference in PSI values between groups to be considered significant.
13. Reference group for differential splicing analysis.
14. Alternative group for differential splicing analysis.
15. Minimum number of reads required to calculate PSI values.
16. True if you want to print PSI values for each sample in the output file.
17. True if you want to perform t-test for differential splicing analysis.
18. True if you want to generate a file of splicing analysis results in excel format.

### 2. Run

``` bash
shiba.py -p 4 config.yaml
```

You are going to use 4 threads for parallelization. You can change the number of threads by changing the `-p` option.

!!! bug "Did you encounter any problems?"

	You can run **Shiba** with the `--verbose` option to see the debug log. This will help you to find the problem.
	```bash
	shiba.py --verbose -p 4 config.yaml
	```
	If you continue to encounter issues, please don't hesitate to [open an issue](https://github.com/Sika-Zheng-Lab/Shiba/issues) on GitHub. The community and developers are here to help!

---

## MameShiba

**MameShiba** is a lightweight version of **Shiba** that can be run on a local machine without Docker or Singularity. It is designed for users who want to perform splicing analysis only and do not need the full functionality of **Shiba**.

### 1. Prepare inputs

`experiment.tsv`: A **tab-separated** text file of sample ID, path to bam files, and groups for differential analysis. This is the same as the input for **Shiba**.

`config.yaml`: A yaml file of the configuration. This is the same as the configuration for **Shiba**.

### 2. Run

Make sure running with `--mame` option.

``` bash
shiba.py --mame -p 4 config.yaml
```

---

## SnakeShiba

A snakemake-based workflow of **Shiba**. This is useful for running **Shiba** on a cluster. Snakemake automatically parallelizes the jobs and manages the dependencies between them.

### 1. Prepare inputs

`experiment.tsv`: A **tab-separated** text file of sample ID, path to fastq files, and groups for differential analysis. This is the same as the input for **Shiba**.

`config.yaml`: A yaml file of the configuration. This is the same as the configuration for **Shiba** but with the addition of the `container` field and without the `only_psi` and `only_psi_group` fields as they are not supported in **SnakeShiba**.

``` yaml
workdir:
  /path/to/workdir # (1)!
container: # This field is required for SnakeShiba
  docker://naotokubota/shiba:v0.8.0 # (2)!
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf # (3)!
experiment_table:
  /path/to/experiment.tsv # (4)!

# Junction read filtering
minimum_anchor_length:
  6 # (5)!
minimum_intron_length:
  70 # (6)!
maximum_intron_length:
  500000 # (7)!
strand:
  XS # (8)!

# PSI calculation
fdr:
  0.05 # (9)!
delta_psi:
  0.1 # (10)!
reference_group:
  Ref # (11)!
alternative_group:
  Alt # (12)!
minimum_reads:
  10 # (13)!
individual_psi:
  True # (14)!
ttest:
  True # (15)!
excel:
  False # (16)!
```

1. The path to the working directory. This is where the output files will be saved. Please make sure that you have write permission to this directory.
2. The Docker image of **Shiba**.
3. The path to the gene annotation file in GTF format.
4. The path to the `experiment.tsv` file.
5. Junctions having a minimum overlap of this much on both ends are reported.
6. Minimum length of the intron sequence.
7. Maximum length of the intron sequence.
8. Strand specificity of RNA library preparation, where the options XS, use XS tags provided by aligner; RF, first-strand; FR, second-strand.
9. Significance threshold for differential splicing analysis.
10. Minimum difference in PSI values between groups to be considered significant.
11. Reference group for differential splicing analysis.
12. Alternative group for differential splicing analysis.
13. Minimum number of reads required to calculate PSI values.
14. True if you want to print PSI values for each sample in the output file.
15. True if you want to perform t-test for differential splicing analysis.
16. True if you want to generate a file of splicing analysis results in excel format.

### 2. Run

Please make sure that you have installed Snakemake and Singularity and cloned the **Shiba** repository on your system.

``` bash
snakemake -s /path/to/Shiba/snakeshiba.smk \
--configfile config.yaml \
--cores 16 \
--use-singularity \
--rerun-incomplete
```
