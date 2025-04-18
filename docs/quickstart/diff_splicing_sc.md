# Differential RNA splicing analysis with single-cell/nucleus RNA-seq data

!!! note "Please install Shiba first"

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

`experiment.tsv`: A **tab-separated** text file of barcode file and STAR solo raw output directory.

``` text
barcode SJ
/path/to/barcodes_run1.tsv /path/to/run1/Solo.out/SJ/raw
/path/to/barcodes_run2.tsv /path/to/run2/Solo.out/SJ/raw
/path/to/barcodes_run3.tsv /path/to/run3/Solo.out/SJ/raw
/path/to/barcodes_run4.tsv /path/to/run4/Solo.out/SJ/raw
```

`barcodes.tsv` is a **tab-separated** text file of barcode and group name like this:

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

``` yaml
workdir:
  /path/to/workdir # (1)!
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf # (2)!
experiment_table:
  /path/to/experiment.tsv # (3)!

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

Docker:

``` bash
cp experiment.tsv config.yaml /path/to/workdir
cd /path/to/workdir
docker run --rm -v $(pwd):$(pwd) -u $(id -u):$(id -g) naotokubota/shiba:v0.5.5 \
python /opt_shiba/Shiba/scshiba.py -p 4 /path/to/workdir/config.yaml
```

Singularity:

``` bash
cp experiment.tsv config.yaml /path/to/workdir
singularity exec docker://naotokubota/shiba:v0.5.5 \
python /opt_shiba/Shiba/scshiba.py -p 4 /path/to/workdir/config.yaml
```

You are going to use 4 threads for parallelization. You can change the number of threads by changing the `-p` option.

!!! note "Binding paths in Singularity"

	When you use Singularity, you do not need to bind any paths as it automatically binds some paths in the host system to the container. In the default configuration, the system default bind points are `$HOME`, `/sys:/sys`, `/proc:/proc`, `/tmp:/tmp`, `/var/tmp:/var/tmp`, `/etc/resolv.conf:/etc/resolv.conf`, `/etc/passwd:/etc/passwd`, and `$PWD`. If files needed to be accessed are not in these paths, you can use the `--bind` option to bind the files to the container.

!!! bug "Did you encounter any problems?"

	You can run **scShiba** with the `--verbose` option to see the debug log. This will help you to find the problem.
	```bash
	docker run --rm -v $(pwd):$(pwd) -u $(id -u):$(id -g) naotokubota/shiba:v0.5.5 \
	python /opt_shiba/Shiba/scshiba.py --verbose -p 4 /path/to/workdir/config.yaml
	```
	If you continue to encounter issues, please don't hesitate to [open an issue](https://github.com/Sika-Zheng-Lab/Shiba/issues) on GitHub. The community and developers are here to help!

---

## SnakeScShiba

A snakemake-based workflow of **scShiba**. This is useful for running **scShiba** on a cluster. Snakemake automatically parallelizes the jobs and manages the dependencies between them.

### 1. Prepare inputs

`experiment.tsv`: A **tab-separated** text file of sample ID, path to fastq files, and groups for differential analysis. This is the same as the input for **scShiba**.

`config.yaml`: A yaml file of the configuration. This is the same as the input for **scShiba** but with the addition of the `container` field.

``` yaml
workdir:
  /path/to/workdir # (1)!
container: # This field is required for SnakeScShiba
  docker://naotokubota/shiba:v0.5.5 # (2)!
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf # (3)!
experiment_table:
  /path/to/experiment.tsv # (4)!

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

Please make sure that you have installed Snakemake and Singularity and cloned the Shiba repository on your system.

``` bash
snakemake -s /path/to/Shiba/snakescshiba.smk \
--configfile config.yaml \
--cores 16 \
--use-singularity \
--rerun-incomplete
```
