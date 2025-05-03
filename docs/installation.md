# Installation

You don't have to build a local environment for running Shiba! We provide a Docker image of **Shiba**, which includes all the dependencies and the **Shiba** software itself. You can run them using [Docker](https://docs.docker.com/get-docker/)/[Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) or [Snakemake](https://snakemake.readthedocs.io/en/stable/).

---

## Docker

``` bash
# Pull the latest image
docker pull naotokubota/shiba:v0.6.0

# Login to the container
docker run -it --rm naotokubota/shiba:v0.6.0 bash

# Run Shiba, for example, to see the help message
docker run -it --rm naotokubota/shiba:v0.6.0 python /opt_shiba/Shiba/shiba.py -h
```

!!! Warning "Memory allocation"

	You may need to allocate more memory to the container if you are using a large dataset. You can do this in the Docker Desktop settings:

	- Go to Docker Desktop settings
	- Click on the "Resources" tab
	- Increase the memory limit as needed
	- Click "Apply & Restart" to save the changes

	![Docker Memory Setting](https://github.com/Sika-Zheng-Lab/Shiba/blob/develop/img/docker_memory_setting.png?raw=true){ align=center width=100% }


---

## Singularity

``` bash
# Pull the latest image
singularity pull docker://naotokubota/shiba:v0.6.0

# Login to the container
singularity shell shiba_v0.6.0.sif

# Run Shiba, for example, to see the help message
singularity exec shiba_v0.6.0.sif python /opt_shiba/Shiba/shiba.py -h
```

---

## Snakemake

You need to install [Snakemake](https://snakemake.readthedocs.io/en/stable/) and clone the **Shiba** GitHub repository on your system:

``` bash
# Clone the Shiba repository
git clone https://github.com/Sika-Zheng-Lab/Shiba.git
```

And please make sure you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) installed on your system as the snakemake workflow uses Singularity to run each step of the pipeline.

---

## Conda (For users who want to run only splicing analysis and do not need full pipeline)

You can also install minimal dependencies for **MameShiba**, a lightweight version of **Shiba** that can be run on a local machine without Docker or Singularity. If you want to perform only splicing analysis, this could be a good option.

The following command will create a conda environment named `mameshiba` with all dependencies installed.

``` bash
conda create -n mameshiba -c conda-forge -c bioconda mameshiba
```
