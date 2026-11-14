# Installation

## Bioconda package

With the `conda-forge` and `bioconda` channels configured:

```bash
conda install shiba
shiba.py --help
scshiba.py --help
```

Run these commands in the activated environment where Shiba is installed.
`shiba.py` and `scshiba.py` are on `PATH`; neither a source checkout nor a
`python ./` prefix is needed. For splicing-only analysis, use
`shiba.py --mame -p 4 config.yaml`. The separate `mameshiba` package also provides
`shiba.py`; there is no `mameshiba` command.

This guide describes the v1.0.0 implementation. Check the version printed by
`shiba.py --help`. If Bioconda still supplies v0.8.2, use the matching
[source checkout](#source-checkout) for the new features. Installing new
dependencies alone does not upgrade Shiba's code. Missing dependencies can be
added with:

```bash
conda install -c conda-forge -c bioconda tosa=1.0.0 \
  r-glmmtmb=1.1.13 r-tmb=1.9.19 r-jsonlite=2.0 openpyxl=3.1
```

## Locating scripts and workflows

Only the two pipeline launchers are installed on `PATH`. For examples below
that use internal scripts, bundled templates, or Snakefiles, set `SHIBA_HOME`
to the directory containing the installed code:

```bash
SHIBA_HOME="$(python -c 'import pathlib, shutil; print(pathlib.Path(shutil.which("shiba.py")).resolve().parent)')"
```

For a source checkout, set `SHIBA_HOME=/absolute/path/to/Shiba` instead.
Run its launchers as `python "$SHIBA_HOME/shiba.py"` or
`python "$SHIBA_HOME/scshiba.py"` in place of `shiba.py` or `scshiba.py` in
this guide. Creating an environment from `environment.yml` installs dependencies,
but does not install these launchers on `PATH`.

For example, internal tools are invoked with
`python "$SHIBA_HOME/src/psi.py" --help`, and workflows with
`snakemake -s "$SHIBA_HOME/snakeshiba.smk" --configfile config.yaml --cores 4`.
Use absolute paths in config and input tables to avoid depending on the directory
from which you launch the command.

## Dependencies for each analysis

The new event-centric factorial regression uses **R/glmmTMB** for model fitting
and **jsonlite** to exchange settings with Python. The table below distinguishes
this requirement from the existing analysis modes.

| Analysis | Additional requirements |
| --- | --- |
| Fisher splicing analysis (optional Welch t-test) | Existing Python dependencies; glmmTMB/jsonlite are not required |
| `stat_method: beta_binomial` | R with glmmTMB and jsonlite |
| Gene-expression differential analysis | Existing R/DESeq2 dependencies, independently of the splicing mode |
| Excel output | `openpyxl` for all workbooks |

The source checkout must contain `src/beta_binomial_glm.R` and the new
`--stat-method` option in `src/psi.py`. Adding R packages does not update an older
Shiba executable. The v0.8.2 package and image contain the earlier counting
pipeline and do not include the factorial regression backend.

## Source checkout

Use the checkout containing the implementation you intend to run. If you are
already working in that checkout, do not clone another copy. Otherwise, clone the
repository and select the branch or release containing the new backend:

```bash
git clone https://github.com/Sika-Zheng-Lab/Shiba.git
cd Shiba
git switch develop # for a released build, select its release tag instead
```

Before continuing, confirm that the selected revision contains the backend:

```bash
test -f src/beta_binomial_glm.R
```

The source-installation commands below are run from this checkout. This explicitly
selects the new source rather than the `shiba.py` executable installed by a
previous Conda package. The new mode is configured as described in the
[beta-binomial regression guide](usage/beta_binomial_regression.md).

## Conda environment for a source checkout

### New Shiba environment

Use the checked-in environment shared by the source install, both Dockerfiles,
and release validation. It includes Excel output and Snakemake:

```bash
conda env create -f environment.yml
conda activate shiba
export PYTHONNOUSERSITE=1
python scripts/verify_runtime.py
```

The regression packages are available as
[`r-glmmtmb`](https://anaconda.org/conda-forge/r-glmmtmb) and
[`r-jsonlite`](https://anaconda.org/conda-forge/r-jsonlite) on conda-forge. Installing
them together lets Conda resolve their R and compiled-library dependencies.
Python 3.12 was used for the regression integration tests; this is a tested
configuration, not a claim that other Python versions cannot work.

Tosa 1.0.0 is used by this checkout for junction and intron-retention boundary
counting in the bulk and single-cell pipelines. The commands above install it
from [Bioconda](https://anaconda.org/bioconda/tosa), making it available on `PATH`
when the Conda environment is active. Tosa is not needed when running `src/psi.py` directly
on existing compatible junction counts and event files.

For optional Excel output:

```bash
python -m pip install openpyxl
```

Run the checked-out pipeline after preparing `config.yaml`:

```bash
python ./shiba.py -p 4 config.yaml
```

### Existing Shiba environment

For reproducible release validation, create a separate environment from this
checkout rather than updating an environment used for older analyses:

```bash
conda env create -n shiba-candidate -f environment.yml
conda activate shiba-candidate
python scripts/verify_runtime.py
```

Keep the environment used for older results available. Adding R packages alone
does not replace an older Shiba executable; run the checked-out scripts.

### MameShiba

Use the same source environment and select splicing-only execution:

```bash
python ./shiba.py --mame -p 4 config.yaml
```

The regression backend is required when `stat_method: beta_binomial` is selected.
`--mame` skips expression, PCA, plotting, and Excel output. The separately
packaged MameShiba dependency set omits the expression backend; the v0.8.2
package does not supply the new factorial backend.

## Installing the R backend with CRAN

If you use a system or module-provided R rather than Conda, install both packages
with the **same Rscript executable** that Shiba will use:

```bash
Rscript --vanilla -e 'install.packages(c("glmmTMB", "jsonlite"), repos="https://cloud.r-project.org")'
```

For a writable, version-specific user library, you can instead use:

```bash
export R_LIBS_USER="$HOME/.local/share/shiba/R-$(Rscript --vanilla -e 'cat(paste(R.version$major, strsplit(R.version$minor, ".", fixed=TRUE)[[1]][1], sep="."))')"
mkdir -p "$R_LIBS_USER"
Rscript --vanilla -e 'install.packages(c("glmmTMB", "jsonlite"), lib=Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")'
```

Keep `R_LIBS_USER` set when running Shiba. A source installation may compile
large C++ dependencies and requires appropriate build tools. Keep glmmTMB and
its TMB/Matrix dependencies compatible; see the
[official glmmTMB installation guide](https://glmmtmb.github.io/glmmTMB/#installation)
for compilation requirements and binary-version mismatch troubleshooting.

By default, Shiba calls `Rscript` from `PATH`. To select another installation,
use `--rscript /path/to/Rscript` with `src/psi.py`, or set this in the pipeline
configuration:

```yaml
rscript: /path/to/Rscript
```

The backend uses `Rscript --vanilla`, so it does not read `.Rprofile`. Install
packages in a library visible to that executable through its standard library
paths or `R_LIBS_USER`, rather than relying on `.Rprofile` to modify `.libPaths()`.
Patsy and Formulaic are not used by the new statistical engine.

## Verify the installation

For the source checkout, run these checks from its root directory in the
environment you will use for analysis. Package users can run `shiba.py --help`
and `python "$SHIBA_HOME/src/psi.py" --help` using the location above. If you configured
a custom `rscript` path, use that executable in the R checks too.

```bash
command -v python
command -v Rscript
python ./shiba.py --help
python ./src/psi.py --help
Rscript --vanilla -e 'stopifnot(requireNamespace("glmmTMB", quietly=TRUE), requireNamespace("jsonlite", quietly=TRUE)); cat(R.version.string, "\n"); print(packageVersion("glmmTMB")); print(packageVersion("jsonlite"))'
```

`src/psi.py --help` should list `--stat-method`, `--formula`, `--coef`, and
`--contrast-file`. For BAM/CRAM counting, also check:

```bash
tosa --version
```

From the source checkout (tests are not included in the release package), run
a fitting check that uses bundled simulated data rather than your BAM files:

```bash
python -m unittest test.test_factorial.TestGlmmTMB.test_known_interaction_and_consistency -v
```

The test must report `ok`, not `skipped`; a skipped test means the R backend was
not available to the test process. Full integration testing used R 4.5.2 with
glmmTMB 1.1.13; selected numerical, event-type and prediction tests also passed
with glmmTMB 1.1.15.2. These are tested versions, not hard minimum versions.
Each analysis records its actual R/glmmTMB versions in `design.json`.

## Docker

Build from the source checkout so that the image contains the new scripts:

```bash
docker build -f docker/Dockerfile -t shiba:factorial .
docker run --rm shiba:factorial tosa --version
docker run --rm shiba:factorial Rscript --vanilla -e 'stopifnot(requireNamespace("glmmTMB", quietly=TRUE), requireNamespace("jsonlite", quietly=TRUE))'
docker run --rm shiba:factorial python /opt/Shiba/src/psi.py --help
```

Both `docker/Dockerfile` and `docker/Dockerfile_develop` install `environment.yml`
with micromamba and copy the checked-out source into `/opt/Shiba`. The release
image targets `linux/amd64`. CI records exact Conda packages, pip packages,
R versions, and the source SHA, then runs every test with no skips. Optional
`shiba2sashimi` is installed separately when needed; it is not required by these
pipelines.
`shiba:factorial` above is a local tag, not a published Docker Hub image.

Mount your working directory and start a shell:

```bash
docker run --rm -it -v "$PWD:/work" -w /work shiba:factorial bash
```

Inside the container, run `shiba.py -p 4 config.yaml`. Input, metadata, contrast,
and output paths in the configuration must be accessible inside the container.
The R packages must be installed in the image; installing them only on the host
does not make them available inside a container.

!!! Warning "Memory allocation"

    Large datasets may require a higher memory limit in Docker Desktop under
    **Settings → Resources**. Apply the new limit before restarting the container.

## Snakemake

Install Snakemake in the environment used to launch the workflow:

```bash
conda install -c conda-forge -c bioconda snakemake
```

Set [`SHIBA_HOME`](#locating-scripts-and-workflows) to the installed code or
source checkout. To run with the activated Conda environment and its dependencies,
omit container execution flags:

```bash
snakemake -s "$SHIBA_HOME/snakeshiba.smk" --configfile config.yaml --cores 4
```

For container execution, install [Apptainer](https://apptainer.org/docs/user/latest/)
on the host and use an image built from the matching source checkout. For
example, convert the Docker image built above into a local SIF image:

```bash
docker save shiba:factorial -o shiba-factorial.tar
apptainer build shiba-factorial.sif docker-archive:shiba-factorial.tar
```

Apptainer supports building SIF images from Docker archives; see its
[container build guide](https://apptainer.org/docs/user/latest/build_a_container.html).
Set the absolute path to the resulting SIF file in `config.yaml`:

```yaml
container: /absolute/path/to/shiba-factorial.sif
```

Then run:

```bash
snakemake -s "$SHIBA_HOME/snakeshiba.smk" --configfile config.yaml --cores 4 --use-apptainer
```

`--use-singularity` is also supported as an alias in the
[Snakemake CLI](https://snakemake.readthedocs.io/en/stable/executing/cli.html#apptainer-singularity).
Ensure that the source checkout and all configured inputs/outputs are mounted
inside the container; use `--apptainer-args` to add site-specific bind mounts.
A locally built Docker tag is not automatically available as a `docker://` registry
image. Replace the example configuration's container URI with your local SIF or
a registry image that contains the new backend.

These runtime setup steps also apply to SnakeScShiba, but the new factorial
regression mode currently targets bulk biological replicate counts in SnakeShiba.

## Fonts for publication figures

Arial is preferred for figures. If it is unavailable, static exports use
Matplotlib's default DejaVu Sans; HTML uses `Arial, sans-serif` and the browser's
fallback font. Optional `SHIBA_FONT_DIR` can point to a directory containing
licensed Arial font files. Fonts are not bundled with this repository.
