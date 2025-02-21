[![GitHub License](https://img.shields.io/github/license/NaotoKubota/Shiba)](https://github.com/NaotoKubota/Shiba/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/801784656.svg)](https://zenodo.org/doi/10.5281/zenodo.11214807)
[![GitHub Release](https://img.shields.io/github/v/release/NaotoKubota/Shiba?style=flat)](https://github.com/NaotoKubota/Shiba/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/NaotoKubota/Shiba)](https://github.com/NaotoKubota/Shiba/releases)
[![Create Release and Build Docker Image](https://github.com/NaotoKubota/Shiba/actions/workflows/release-docker-build-push.yaml/badge.svg)](https://github.com/NaotoKubota/Shiba/actions/workflows/release-docker-build-push.yaml)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/shiba)](https://hub.docker.com/r/naotokubota/shiba)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/shiba)](https://hub.docker.com/r/naotokubota/shiba)

# Shiba (v0.5.2) <img src="img/Shiba_icon.png" width=40% align="right">

A versatile computational method for robust identification of differential RNA splicing. Shiba/scShiba can quantify and identify differential splicing events (DSEs) from short-read bulk RNA-seq data and single-cell RNA-seq data. Shiba and scShiba are also implemented as [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows, SnakeShiba and SnakeScShiba, respectively.

See [CHANGELOG.md](https://github.com/NaotoKubota/Shiba/blob/main/CHANGELOG.md) for the latest updates.

## Overview

Shiba comprises four main steps:
1. **Transcript assembly**: Assemble transcripts from RNA-seq reads using [StringTie2](https://github.com/skovaka/stringtie2)
2. **Splicing event identification**: Identify alternative mRNA splicing events from assembled transcripts
3. **Read counting**: Count reads mapped to each splicing event using [RegTools](https://github.com/griffithlab/regtools) and [featureCounts](https://subread.sourceforge.net/)
4. **Statistical analysis**: Identify DSEs based on Fisher's exact test

<img src="img/Shiba_overview.png" width=75%>

## Installation

A docker image is available at [Docker Hub](https://hub.docker.com/r/naotokubota/shiba).

```bash
docker pull naotokubota/shiba:v0.5.2
```

## Usage

Manual for Shiba is available at [https://sika-zheng-lab.github.io/Shiba/](https://sika-zheng-lab.github.io/Shiba/).

***Shiba***

```bash
python shiba.py -p 32 config.yaml
```

***SnakeShiba***, Snakemake-based workflow of Shiba

```bash
snakemake -s snakeshiba.smk --configfile config.yaml --cores 32 --use-singularity
```

***scShiba***, a single-cell RNA-seq version of Shiba

```bash
python scshiba.py -p 32 config.yaml
```

***SnakeScShiba***, Snakemake-based workflow of scShiba

```bash
snakemake -s snakescshiba.smk --configfile config.yaml --cores 32 --use-singularity
```

## Citation

Kubota N, Chen L, Zheng S. (2024). [Shiba: A unified computational method for robust identification of differential RNA splicing across platforms](https://www.biorxiv.org/content/10.1101/2024.05.30.596331v1). *bioRxiv* 2024.05.30.596331

## Authors

- Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
- Liang Chen ([0000-0001-6164-4553](https://orcid.org/0000-0001-6164-4553))
- Sika Zheng ([0000-0002-0573-4981](https://orcid.org/0000-0002-0573-4981))
