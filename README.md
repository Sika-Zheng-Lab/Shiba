[![GitHub License](https://img.shields.io/github/license/Sika-Zheng-Lab/Shiba)](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/801784656.svg)](https://zenodo.org/doi/10.5281/zenodo.11214807)
[![GitHub Release](https://img.shields.io/github/v/release/Sika-Zheng-Lab/Shiba?style=flat)](https://github.com/Sika-Zheng-Lab/Shiba/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/Sika-Zheng-Lab/Shiba)](https://github.com/Sika-Zheng-Lab/Shiba/releases)
[![Create Release and Build Docker Image](https://github.com/Sika-Zheng-Lab/Shiba/actions/workflows/release-docker-build-push.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/Shiba/actions/workflows/release-docker-build-push.yaml)
[![Conda](https://img.shields.io/conda/v/bioconda/mameshiba?color=3EB049)](https://anaconda.org/bioconda/mameshiba)
[![Conda - Downloads](https://img.shields.io/conda/dn/bioconda/mameshiba?label=Conda%20-%20Downloads&color=3EB049)](https://anaconda.org/bioconda/mameshiba)
[![Docker](https://img.shields.io/docker/v/naotokubota/shiba?color=blue&label=Docker)](https://hub.docker.com/r/naotokubota/shiba)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/shiba)](https://hub.docker.com/r/naotokubota/shiba)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/shiba)](https://hub.docker.com/r/naotokubota/shiba)
[![NAR](https://img.shields.io/badge/NAR-10.1093/nar/gkaf098-0B3B58)](https://academic.oup.com/nar/article/53/4/gkaf098/8042001)

# Shiba (v0.6.0) <img src="https://raw.githubusercontent.com/Sika-Zheng-Lab/Shiba/main/img/Shiba_icon.png" width=40% align="right">

A versatile computational method for systematic identification of differential RNA splicing. Shiba/scShiba can quantify and identify differential splicing events (DSEs) from bulk RNA-seq data and single-cell RNA-seq data. Shiba and scShiba are also implemented as [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows, SnakeShiba and SnakeScShiba, respectively.

See [CHANGELOG.md](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/CHANGELOG.md) for the latest updates.

> [!IMPORTANT]
> Shiba now supports **long-read RNA-seq** data analysis! Please check the [manual](https://sika-zheng-lab.github.io/Shiba/quickstart/diff_splicing_bulk/#shiba) for more details.

## Overview

Shiba comprises four main steps:
1. **Transcript assembly**: Assemble transcripts from RNA-seq reads using [StringTie2](https://github.com/skovaka/stringtie2)
2. **Splicing event identification**: Identify alternative mRNA splicing events from assembled transcripts
3. **Read counting**: Count reads mapped to each splicing event using [RegTools](https://github.com/griffithlab/regtools) and [featureCounts](https://subread.sourceforge.net/)
4. **Statistical analysis**: Identify DSEs based on Fisher's exact test

<img src="https://raw.githubusercontent.com/Sika-Zheng-Lab/Shiba/main/img/Shiba_overview.png" width=75%>

## Installation

### Conda

```bash
conda create -n shiba -c conda-forge -c bioconda shiba
```

If you want to perform only splicing analysis, you can install minimal dependencies and run **MameShiba**, a lightweight version of Shiba.

```bash
conda create -n mameshiba -c conda-forge -c bioconda mameshiba
```

### Docker

```bash
docker pull naotokubota/shiba:v0.6.0
```

## Usage

Manual for Shiba is available at [https://sika-zheng-lab.github.io/Shiba/](https://sika-zheng-lab.github.io/Shiba/).

***Shiba***

```bash
shiba.py -p 4 config.yaml
```

***MameShiba***, a lightweight version of Shiba

```bash
shiba.py --mame -p 4 config.yaml
```

***SnakeShiba***, Snakemake-based workflow of Shiba

```bash
snakemake -s snakeshiba.smk --configfile config.yaml --cores 4 --use-singularity
```

***scShiba***, a single-cell RNA-seq version of Shiba

```bash
scshiba.py -p 4 config.yaml
```

***SnakeScShiba***, Snakemake-based workflow of scShiba

```bash
snakemake -s snakescshiba.smk --configfile config.yaml --cores 4 --use-singularity
```

## Visualization

Do you want to visualize the results of Shiba analysis? Try 🐕 [shiba2sashimi](https://github.com/Sika-Zheng-Lab/shiba2sashimi) 🍣 !

<img src="https://raw.githubusercontent.com/Sika-Zheng-Lab/shiba2sashimi/main/img/sashimi_example.png" width=80%>

## Contributing

Thank you for wanting to improve Shiba! If you have any bugs or questions, feel free to [open an issue](https://github.com/Sika-Zheng-Lab/Shiba/issues) or pull request.

## Citation

Kubota N, Chen L, Zheng S. [Shiba: a versatile computational method for systematic identification of differential RNA splicing across platforms](https://academic.oup.com/nar/article/53/4/gkaf098/8042001). *Nucleic Acids Research*  53(4), 2025, gkaf098.

## Authors

- Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
- Liang Chen ([0000-0001-6164-4553](https://orcid.org/0000-0001-6164-4553))
- Sika Zheng ([0000-0002-0573-4981](https://orcid.org/0000-0002-0573-4981))
