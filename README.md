[![GitHub License](https://img.shields.io/github/license/Sika-Zheng-Lab/Shiba)](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/801784656.svg)](https://zenodo.org/doi/10.5281/zenodo.11214807)
[![GitHub Release](https://img.shields.io/github/v/release/Sika-Zheng-Lab/Shiba?style=flat)](https://github.com/Sika-Zheng-Lab/Shiba/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/Sika-Zheng-Lab/Shiba)](https://github.com/Sika-Zheng-Lab/Shiba/releases)
[![Create Release and Build Docker Image](https://github.com/Sika-Zheng-Lab/Shiba/actions/workflows/release-docker-build-push.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/Shiba/actions/workflows/release-docker-build-push.yaml)
[![Tests](https://github.com/Sika-Zheng-Lab/Shiba/actions/workflows/test.yaml/badge.svg)](https://github.com/Sika-Zheng-Lab/Shiba/actions/workflows/test.yaml)
[![codecov](https://codecov.io/gh/Sika-Zheng-Lab/Shiba/branch/develop/graph/badge.svg?token=GCBM8OWSLY)](https://codecov.io/gh/Sika-Zheng-Lab/Shiba)
[![Conda](https://img.shields.io/conda/v/bioconda/shiba?color=3EB049)](https://anaconda.org/bioconda/shiba)
[![Conda - Shiba - Downloads](https://img.shields.io/conda/dn/bioconda/shiba?label=Conda%20-%20Shiba%20-%20Downloads&color=3EB049)](https://anaconda.org/bioconda/shiba)
[![Conda - MameShiba - Downloads](https://img.shields.io/conda/dn/bioconda/mameshiba?label=Conda%20-%20MameShiba%20-%20Downloads&color=3EB049)](https://anaconda.org/bioconda/mameshiba)
[![Docker](https://img.shields.io/docker/v/naotokubota/shiba?color=blue&label=Docker)](https://hub.docker.com/r/naotokubota/shiba)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/shiba)](https://hub.docker.com/r/naotokubota/shiba)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/shiba)](https://hub.docker.com/r/naotokubota/shiba)
[![NAR](https://img.shields.io/badge/NAR-10.1093/nar/gkaf098-0B3B58)](https://academic.oup.com/nar/article/53/4/gkaf098/8042001)

# Shiba (v0.8.2) <img src="https://raw.githubusercontent.com/Sika-Zheng-Lab/Shiba/main/img/Shiba_icon.png" width=40% align="right">

A versatile computational method for systematic identification of differential RNA splicing. Shiba/scShiba can quantify and identify differential splicing events (DSEs) from bulk RNA-seq data and single-cell RNA-seq data. Shiba and scShiba are also implemented as [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows, SnakeShiba and SnakeScShiba, respectively.

The scShiba and SnakeScShiba code in this checkout counts junction and intron-retention boundary reads directly from tagged BAM/CRAM files with Tosa. See the [single-cell input guide](docs/quickstart/diff_splicing_sc.md) for the new experiment table.

See [CHANGELOG.md](https://github.com/Sika-Zheng-Lab/Shiba/blob/main/CHANGELOG.md) for the latest updates.

> [!IMPORTANT]
> Shiba now supports **long-read RNA-seq** data analysis! Please check the [manual](https://sika-zheng-lab.github.io/Shiba/quickstart/diff_splicing_bulk/#shiba) for more details.

## Overview

Shiba comprises four main steps:
1. **Transcript assembly**: Assemble transcripts from RNA-seq reads using [StringTie](https://github.com/gpertea/stringtie)
2. **Splicing event identification**: Identify alternative mRNA splicing events from assembled transcripts
3. **Read counting**: Count splice junction and intron retention boundary reads using [Tosa](https://github.com/Sika-Zheng-Lab/Tosa)
4. **Statistical analysis**: Identify DSEs using Fisher's exact test or beta-binomial regression

**Beta-binomial analysis**: Shiba and SnakeShiba support optional event-centric
beta-binomial regression (`stat_method: beta_binomial`) for factorial experimental
designs. This mode models junction counts across biological replicates with
R/glmmTMB, accounting for overdispersion and supporting covariates, interaction
terms (for example, genotype × treatment), and user-defined contrasts. It reports
FDR-adjusted event-level results and model-based PSI predictions and effect sizes.
See the [beta-binomial regression guide](docs/usage/beta_binomial_regression.md)
for configuration and examples. Set `stat_method: beta_binomial` plus
`formula`, metadata, and `coef` or `contrast_file` in the configuration.
Set `stat_method: fisher` to use Fisher testing instead. The old supplementary
`beta_binomial` key and `psi.py -b/--beta-binomial` have been removed; delete
that key from existing configurations.

<img src="https://raw.githubusercontent.com/Sika-Zheng-Lab/Shiba/main/img/Shiba_overview.png" width=75%>

## Installation

Install Shiba from Bioconda (with the `conda-forge` and `bioconda` channels configured):

```bash
conda install shiba
```

If the new dependencies are missing from an existing environment, install them with:

```bash
conda install -c conda-forge -c bioconda tosa=1.0.0 \
  r-glmmtmb=1.1.13 r-tmb=1.9.19 r-jsonlite=2.0 openpyxl=3.1
```

If Bioconda has not yet published Shiba v1.0.0, these dependencies alone do not
upgrade the Shiba code. Use the [source installation instructions](docs/installation.md#source-checkout)
with the v1.0.0 release checkout and its `environment.yml` until the package is available.
The source environment also supplies the matching Python, plotting, and expression-analysis dependencies.

See the [installation guide](docs/installation.md) for other installation options.

## Usage

See the documentation for [pipeline commands](docs/quickstart/running.md),
[bulk RNA-seq](docs/quickstart/diff_splicing_bulk.md),
[single-cell RNA-seq](docs/quickstart/diff_splicing_sc.md), and
[factorial beta-binomial analysis with example inputs](docs/quickstart/beta_binomial.md).

## Visualization

Do you want to visualize the results of Shiba analysis? Try 🐕 [shiba2sashimi](https://github.com/Sika-Zheng-Lab/shiba2sashimi) 🍣 !

<img src="https://raw.githubusercontent.com/Sika-Zheng-Lab/shiba2sashimi/main/img/sashimi_example.png" width=80%>

## Contributing

Thank you for wanting to improve Shiba! If you have any bugs or questions, feel free to [open an issue](https://github.com/Sika-Zheng-Lab/Shiba/issues) or pull request.

### Preview the documentation

The user guide is built with [Zensical](https://zensical.org/), using the existing
`mkdocs.yml` configuration and the classic theme. From the repository root, use
Python 3.12 to create a separate documentation environment:

```bash
python3.12 -m venv .venv-docs
source .venv-docs/bin/activate
python -m pip install -r requirements-docs.txt
zensical serve
```

Open <http://localhost:8000/Shiba/>. Changes to `docs/` are rebuilt automatically;
press Ctrl+C to stop the preview. To run the same build and link validation as CI:

```bash
zensical build --clean --strict
```

The generated site is written to `site/`. If running on a remote server, forward
port 8000 to your computer (for example, `ssh -L 8000:127.0.0.1:8000 user@host`)
and open the same local URL. If port 8000 is busy, use
`zensical serve --dev-addr 127.0.0.1:8001` and forward/open port 8001 instead.

Release CI publishes the verified build to the existing `gh-pages` branch with
`ghp-import` after the release succeeds. Local preview and build commands do not
publish the site.

## Citation

Kubota N, Chen L, Zheng S. [Shiba: a versatile computational method for systematic identification of differential RNA splicing across platforms](https://academic.oup.com/nar/article/53/4/gkaf098/8042001). *Nucleic Acids Research*  53(4), 2025, gkaf098.

## Authors

- Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
- Liang Chen ([0000-0001-6164-4553](https://orcid.org/0000-0001-6164-4553))
- Sika Zheng ([0000-0002-0573-4981](https://orcid.org/0000-0002-0573-4981))
