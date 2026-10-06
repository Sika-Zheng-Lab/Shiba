# Installation

The new event-centric factorial regression uses **R/glmmTMB** for model fitting
and **jsonlite** to exchange settings with Python. The table below distinguishes
this requirement from the existing analysis modes.

| Analysis | Additional requirements |
| --- | --- |
| Legacy splicing, including the older `beta_binomial: True` supplementary test | Existing Python dependencies; glmmTMB/jsonlite are not required |
| `stat_method: beta-binomial` | R with glmmTMB and jsonlite |
| Gene-expression differential analysis | Existing R/DESeq2 dependencies, independently of the splicing mode |
| Excel output | `styleframe==4.2` for legacy output; `openpyxl` for the new regression workbook |

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
```

Before continuing, confirm that the selected revision contains the backend:

```bash
test -f src/beta_binomial_glm.R
```

All `python ./...` commands below are run from this checkout. This explicitly
selects the new source rather than the `shiba.py` executable installed by a
previous Conda package. The new mode is configured as described in the
[beta-binomial regression guide](usage/beta_binomial_regression.md).

## Conda

### New Shiba environment

Use the Shiba package to install the base pipeline dependencies, then run the
checked-out source:

```bash
conda create -n shiba -c conda-forge -c bioconda \
  python=3.12 shiba r-glmmtmb r-jsonlite rust
conda activate shiba
cargo install tosa --version 1.0.0 --locked
export PATH="${CARGO_HOME:-$HOME/.cargo}/bin:$PATH"
```

The regression packages are available as
[`r-glmmtmb`](https://anaconda.org/conda-forge/r-glmmtmb) and
[`r-jsonlite`](https://anaconda.org/conda-forge/r-jsonlite) on conda-forge. Installing
them together lets Conda resolve their R and compiled-library dependencies.
Python 3.12 was used for the regression integration tests; this is a tested
configuration, not a claim that other Python versions cannot work.

Tosa 1.0.0 is used by this checkout for junction and intron-retention boundary
counting in the bulk and single-cell pipelines. The command above installs it
with Cargo and adds Cargo's executable directory to `PATH`; retain that PATH
setting in future sessions. Tosa is not needed when running `src/psi.py` directly
on existing compatible junction counts and event files.

For optional Excel output:

```bash
python -m pip install styleframe==4.2 openpyxl
```

Run the checked-out pipeline after preparing `config.yaml`:

```bash
python ./shiba.py -p 4 config.yaml
```

### Existing Shiba environment

Add the R backend to the environment used to run Shiba:

```bash
conda activate shiba
conda install -c conda-forge r-glmmtmb r-jsonlite
```

Use `python ./shiba.py` or `python ./src/psi.py` from the new checkout after the
installation. The additional packages alone do not replace an older installed
Shiba script. If counting from BAM/CRAM files, also install Tosa as shown above.

### MameShiba

For splicing-only analysis, use the MameShiba package for the base dependencies:

```bash
conda create -n mameshiba -c conda-forge -c bioconda \
  python=3.12 mameshiba r-glmmtmb r-jsonlite rust
conda activate mameshiba
cargo install tosa --version 1.0.0 --locked
export PATH="${CARGO_HOME:-$HOME/.cargo}/bin:$PATH"
python ./shiba.py --mame -p 4 config.yaml
```

The new regression backend is required even in MameShiba when
`stat_method: beta-binomial` is selected. The full pipeline's expression, PCA,
and plotting steps are skipped by `--mame`; its existing Excel behavior is
unchanged.

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

Run these checks in the environment you will use for analysis. If you configured
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

For a fitting check that uses bundled simulated data rather than your BAM files:

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
docker run --rm shiba:factorial Rscript --vanilla -e 'stopifnot(requireNamespace("glmmTMB", quietly=TRUE), requireNamespace("jsonlite", quietly=TRUE))'
docker run --rm shiba:factorial python /opt/Shiba/src/psi.py --help
```

Both `docker/Dockerfile` and `docker/Dockerfile_develop` install glmmTMB/jsonlite,
install Tosa 1.0.0, and copy the checked-out Shiba source into `/opt/Shiba`.
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

To run with the activated Conda environment and its dependencies, omit container
execution flags:

```bash
snakemake -s snakeshiba.smk --configfile config.yaml --cores 4
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
snakemake -s snakeshiba.smk --configfile config.yaml --cores 4 --use-apptainer
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
