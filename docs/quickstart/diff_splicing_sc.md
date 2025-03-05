# Differential RNA splicing analysis with single-cell/nucleus RNA-seq data

!!! note

	You need to install a Docker image of **Shiba** (and clone the **Shiba** GitHub repository to run **SnakeScShiba**). If you don't have them installed, please follow the instructions in the [Installation](../installation.md) section.

---

## Before you start

- Perform mapping of sc(sn)RNA-seq reads to the reference genome using [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).
  - You can download a test input file mapped by STARsolo on the mouse genome from [here](https://zenodo.org/records/14976391).
- Download a gene annotataion file of your interest in GTF format.

Here is an example code for downloading a mouse gene annotation file (Ensembl 102):

``` bash
wget https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gzip -d Mus_musculus.GRCm38.102.gtf.gz
```

---

## scShiba

### 1. Prepare inputs

`experiment.tsv`: A tab-separated text file of barcode file and STAR solo raw output directory.

``` text
barcode<tab>SJ
/path/to/barcodes_run1.tsv<tab>/path/to/run1/Solo.out/SJ/raw
/path/to/barcodes_run2.tsv<tab>/path/to/run2/Solo.out/SJ/raw
/path/to/barcodes_run3.tsv<tab>/path/to/run3/Solo.out/SJ/raw
/path/to/barcodes_run4.tsv<tab>/path/to/run4/Solo.out/SJ/raw
```

`barcodes.tsv` is a tab-separated text file of barcode and group name like this:

``` text
barcode<tab>group
TTTGTTGTCCACACCT<tab>Cluster-1
TCAAGACCACTACAGT<tab>Cluster-1
TATTTCGGTACAGTAA<tab>Cluster-1
ATCCTATGTTAATCGC<tab>Cluster-1
ATCGATGAGTTTCTTC<tab>Cluster-2
ATCGATGGTCTTGCTC<tab>Cluster-2
TATGTTCGTCAGGCAA<tab>Cluster-2
ATCGCCTAGACTCGAG<tab>Cluster-2
...
```

Please replace `<tab>` with a tab character.

`config.yaml`: A yaml file of the configuration.

``` yaml
workdir:
  /path/to/workdir
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf
experiment_table:
  /path/to/experiment.tsv

# PSI calculation
only_psi:
  False
fdr:
  0.05
delta_psi:
  0.1
reference_group:
  Cluster-1
alternative_group:
  Cluster-2
minimum_reads:
  10
excel:
  True
```

You can generate a file of splicing analysis results in excel format by setting `excel` to `True`.

### 2. Run

Docker:

``` bash
cp experiment.tsv config.yaml /path/to/workdir
cd /path/to/workdir
docker run --rm -v $(pwd):$(pwd) -u $(id -u):$(id -g) naotokubota/shiba:v0.5.2 \
python /opt_shiba/Shiba/scshiba.py -p 32 /path/to/workdir/config.yaml
```

Singularity:

``` bash
cp experiment.tsv config.yaml /path/to/workdir
singularity exec docker://naotokubota/shiba:v0.5.2 \
python /opt_shiba/Shiba/scshiba.py -p 32 /path/to/workdir/config.yaml
```

!!! note

	When you use Singularity, you do not need to bind any paths as it automatically binds some paths in the host system to the container. In the default configuration, the system default bind points are `$HOME`, `/sys:/sys`, `/proc:/proc`, `/tmp:/tmp`, `/var/tmp:/var/tmp`, `/etc/resolv.conf:/etc/resolv.conf`, `/etc/passwd:/etc/passwd`, and `$PWD`. If files needed to be accessed are not in these paths, you can use the `--bind` option to bind the files to the container.

---

## SnakeScShiba

A snakemake-based workflow of **scShiba**. This is useful for running **scShiba** on a cluster. Snakemake automatically parallelizes the jobs and manages the dependencies between them.

### 1. Prepare inputs

`experiment.tsv`: A tab-separated text file of sample ID, path to fastq files, and groups for differential analysis. This is the same as the input for **scShiba**.

`config.yaml`: A yaml file of the configuration. This is the same as the input for **scShiba** but with the addition of the `container` field.

``` yaml
workdir:
  /path/to/workdir
container: # This field is required for SnakeScShiba
  docker://naotokubota/shiba:v0.5.2
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf
experiment_table:
  /path/to/experiment.tsv

# PSI calculation
only_psi:
  False
fdr:
  0.05
delta_psi:
  0.1
reference_group:
  Cluster-1
alternative_group:
  Cluster-2
minimum_reads:
  10
excel:
  False
```

You can generate a file of splicing analysis results in excel format by setting `excel` to `True`.

### 2. Run

Please make sure that you have installed Snakemake and Singularity and cloned the Shiba repository on your system.

``` bash
snakemake -s /path/to/Shiba/snakescshiba.smk \
--configfile config.yaml \
--cores 16 \
--use-singularity \
--rerun-incomplete
```
