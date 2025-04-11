# Differential RNA splicing analysis with bulk RNA-seq data

!!! note

	You need to install a Docker image of **Shiba** (and clone the **Shiba** GitHub repository to run **SnakeShiba**). If you don't have them installed, please follow the instructions in the [Installation](../installation.md) section.

---

## Before you start

- Perform mapping of RNA-seq reads to the reference genome and generate bam files by software such as [STAR](https://github.com/alexdobin/STAR) and [HISAT2](https://daehwankimlab.github.io/hisat2/).
  - You can download test RNA-seq bam files with their index (two replicates for reference and alternative groups) mapped by STAR on the mouse genome from [here](https://zenodo.org/records/14976391).
- Download a gene annotataion file of your interest in GTF format.

Here is an example code for downloading a mouse gene annotation file (Ensembl 102):

``` bash
wget https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gzip -d Mus_musculus.GRCm38.102.gtf.gz
```

---

## Shiba

### 1. Prepare inputs

`experiment.tsv`: A tab-separated text file of sample ID, path to bam files, and groups for differential analysis.

``` text
sample<tab>bam<tab>group
Ref_1<tab>/path/to/workdir/bam/Ref_1.bam<tab>Ref
Ref_2<tab>/path/to/workdir/bam/Ref_2.bam<tab>Ref
Alt_1<tab>/path/to/workdir/bam/Alt_1.bam<tab>Alt
Alt_2<tab>/path/to/workdir/bam/Alt_2.bam<tab>Alt
```

Please put bam files with their index files (`.bai`) in the `path/to/workdir/bam` directory and replace `<tab>` with a tab character.

`config.yaml`: A yaml file of the configuration.

``` yaml
workdir:
  /path/to/workdir
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf
experiment_table:
  /path/to/experiment.tsv
unannotated:
  True

# Junction read filtering
minimum_anchor_length:
  6
minimum_intron_length:
  70
maximum_intron_length:
  500000
strand:
  XS

# PSI calculation
only_psi:
  False # True if you want to calculate PSI values for each sample only
only_psi_group:
  False # True if you want to calculate PSI values for each group only (Overrides `only_psi`)
fdr:
  0.05
delta_psi:
  0.1
reference_group:
  Ref
alternative_group:
  Alt
minimum_reads:
  10
individual_psi:
  True
ttest:
  True
excel:
  False # True if you want to generate a file of splicing analysis results in excel format
```

### 2. Run

Docker:

``` bash
cp experiment.tsv config.yaml /path/to/workdir
cd /path/to/workdir
docker run --rm -v $(pwd):$(pwd) -u $(id -u):$(id -g) naotokubota/shiba:v0.5.5 \
python /opt_shiba/Shiba/shiba.py -p 4 /path/to/workdir/config.yaml
```

Singularity:

``` bash
cp experiment.tsv config.yaml /path/to/workdir
singularity exec docker://naotokubota/shiba:v0.5.5 \
python /opt_shiba/Shiba/shiba.py -p 4 /path/to/workdir/config.yaml
```

You are going to use 4 threads for parallelization. You can change the number of threads by changing the `-p` option.

!!! note

	When you use Singularity, you do not need to bind any paths as it automatically binds some paths in the host system to the container. In the default configuration, the system default bind points are `$HOME`, `/sys:/sys`, `/proc:/proc`, `/tmp:/tmp`, `/var/tmp:/var/tmp`, `/etc/resolv.conf:/etc/resolv.conf`, `/etc/passwd:/etc/passwd`, and `$PWD`. If files needed to be accessed are not in these paths, you can use the `--bind` option to bind the files to the container.

!!! tip

	Did you encounter any problems?
	You can run **Shiba** with the `--verbose` option to see the debug log. This will help you to find the problem.
	```bash
	docker run --rm -v $(pwd):$(pwd) -u $(id -u):$(id -g) naotokubota/shiba:v0.5.5 \
	python /opt_shiba/Shiba/shiba.py --verbose -p 4 /path/to/workdir/config.yaml
	```
	If you continue to encounter issues, please don't hesitate to [open an issue](https://github.com/Sika-Zheng-Lab/Shiba/issues) on GitHub. The community and developers are here to help!

---

## SnakeShiba

A snakemake-based workflow of **Shiba**. This is useful for running **Shiba** on a cluster. Snakemake automatically parallelizes the jobs and manages the dependencies between them.

### 1. Prepare inputs

`experiment.tsv`: A tab-separated text file of sample ID, path to fastq files, and groups for differential analysis. This is the same as the input for **Shiba**.

`config.yaml`: A yaml file of the configuration. This is the same as the configuration for **Shiba** but with the addition of the `container` field and without the `only_psi` and `only_psi_group` fields as they are not supported in **SnakeShiba**.

``` yaml
workdir:
  /path/to/workdir
container: # This field is required for SnakeShiba
  docker://naotokubota/shiba:v0.5.5
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf
experiment_table:
  /path/to/experiment.tsv
unannotated:
  True

# Junction read filtering
minimum_anchor_length:
  6
minimum_intron_length:
  70
maximum_intron_length:
  500000
strand:
  XS

# PSI calculation
fdr:
  0.05
delta_psi:
  0.1
reference_group:
  Ref
alternative_group:
  Alt
minimum_reads:
  10
individual_psi:
  True
ttest:
  True
excel:
  False # True if you want to generate a file of splicing analysis results in excel format
```

### 2. Run

Please make sure that you have installed Snakemake and Singularity and cloned the **Shiba** repository on your system.

``` bash
snakemake -s /path/to/Shiba/snakeshiba.smk \
--configfile config.yaml \
--cores 16 \
--use-singularity \
--rerun-incomplete
```
