###############################################################
# Dockerfile to build container images for Shiba (v0.6.3)
# Based on continuumio/miniconda3:23.10.0-1
################################################################

FROM continuumio/miniconda3:23.10.0-1

# File Author / Maintainer
LABEL maintainer="Naoto Kubota <naotok@ucr.edu>"

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies first
RUN apt-get -qq update && \
	apt-get -qq -y install \
    build-essential gcc g++ make cmake libtool texinfo dpkg-dev pkg-config \
    libgsl-dev wget locate less vim zlib1g-dev bzip2 lzma curl r-base \
    libboost-dev libcurl4-openssl-dev libboost-all-dev libbz2-dev liblzma-dev \
    libpcre3 libpcre3-dev

# Install mamba
RUN conda install -c conda-forge mamba

# Install Shiba with conda
RUN mamba install -c conda-forge -c bioconda -c defaults mameshiba plotly==5.13.0

# Install R (version 4.1.3)
RUN mamba install -c conda-forge -c bioconda r-base==4.1.3 bioconductor-deseq2=1.34.0 r-locfit==1.5_9.4 r-data.table==1.14.2

# Install python packages
RUN pip install --upgrade pip && \
    pip install styleframe==4.1

# Install shiba2sashimi
RUN pip install shiba2sashimi

# Install Shiba repository
RUN cd /opt && \
	git clone https://github.com/Sika-Zheng-Lab/Shiba.git -b v0.6.3 && \
	ln -sf /opt/Shiba/shiba.py /opt/conda/bin/shiba.py && \
	ln -sf /opt/Shiba/scshiba.py /opt/conda/bin/scshiba.py

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
