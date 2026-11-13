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

The event-centric factorial beta-binomial mode requires **R, glmmTMB, and
jsonlite**, in addition to Shiba's existing dependencies. Use a source checkout
that contains `src/beta_binomial_glm.R`; installing these R packages alone does
not add the new mode to an older Shiba executable. The v0.8.2 package/image uses
the earlier counting pipeline and does not contain this regression backend.

### Conda

Use the environment shared by source installs, containers, and release checks:

```bash
git clone https://github.com/Sika-Zheng-Lab/Shiba.git
cd Shiba
git switch develop # select a release tag instead when installing that release
conda env create -f environment.yml
conda activate shiba
python scripts/verify_runtime.py
python ./shiba.py -p 4 config.yaml
```

`python ./shiba.py` selects this checkout instead of an older executable on
`PATH`. The environment includes Tosa, the R regression backend, gene-expression
dependencies, Snakemake, and Excel support. For splicing-only **MameShiba**, run
`python ./shiba.py --mame -p 4 config.yaml`; regression still requires glmmTMB/jsonlite.
See the [installation guide](docs/installation.md) for existing environments and
optional R setup.

### Docker

Build the same checkout locally (the supported release image platform is
`linux/amd64`):

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t shiba:local .
docker run --rm shiba:local python /opt/Shiba/scripts/verify_runtime.py
```

Published images use the **`naotokubota/shiba`** repository. Select an existing
release tag that matches the scripts and input format; the v0.8.2 image does
not contain the Tosa/factorial workflow described here.

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

## Citation

Kubota N, Chen L, Zheng S. [Shiba: a versatile computational method for systematic identification of differential RNA splicing across platforms](https://academic.oup.com/nar/article/53/4/gkaf098/8042001). *Nucleic Acids Research*  53(4), 2025, gkaf098.

## Authors

- Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
- Liang Chen ([0000-0001-6164-4553](https://orcid.org/0000-0001-6164-4553))
- Sika Zheng ([0000-0002-0573-4981](https://orcid.org/0000-0002-0573-4981))
