# Installation

## Conda

The following command will create a conda environment named `shiba` with all dependencies installed.

``` bash
conda create -n shiba -c conda-forge -c bioconda shiba
conda activate shiba # Activate the conda environment
pip install styleframe==4.1 # optional, for generating outputs in Excel format.
```

You can also install minimal dependencies for **MameShiba**, a lightweight version of **Shiba** . If you want to perform only splicing analysis, this could be a good option. The following command will create a conda environment named `mameshiba` with minimal dependencies installed.

``` bash
conda create -n mameshiba -c conda-forge -c bioconda mameshiba
```

---

## Docker

We provide a Docker image for **Shiba**. You can use the following command to pull the latest image from Docker Hub:

``` bash
# Pull the latest image
docker pull naotokubota/shiba:v0.6.3

# Login to the container
docker run -it --rm naotokubota/shiba:v0.6.3 bash

# Run Shiba, for example, to see the help message
docker run -it --rm naotokubota/shiba:v0.6.3 shiba.py -h
```

!!! Warning "Memory allocation"

	You may need to allocate more memory to the container if you are using a large dataset. You can do this in the Docker Desktop settings:

	- Go to Docker Desktop settings
	- Click on the "Resources" tab
	- Increase the memory limit as needed
	- Click "Apply & Restart" to save the changes

	![Docker Memory Setting](https://github.com/Sika-Zheng-Lab/Shiba/blob/develop/img/docker_memory_setting.png?raw=true){ align=center width=100% }


---

## Snakemake

You need to install [Snakemake](https://snakemake.readthedocs.io/en/stable/) and clone the **Shiba** GitHub repository on your system:

``` bash
# Clone the Shiba repository
git clone https://github.com/Sika-Zheng-Lab/Shiba.git
```

And please make sure you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) installed on your system as the snakemake workflow uses Singularity to run each step of the pipeline.
