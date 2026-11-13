###############################################################
# Dockerfile to build container images for Shiba
# Based on mambaorg/micromamba:2.5
################################################################

FROM mambaorg/micromamba:2.5

# File Author / Maintainer
LABEL maintainer="Naoto Kubota <naotok@ucr.edu>"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH

USER root

# Install dependencies first
RUN apt-get -qq update && \
	apt-get -qq -y install \
    build-essential gcc g++ make cmake libtool texinfo dpkg-dev pkg-config \
    libgsl-dev wget locate less vim zlib1g-dev bzip2 lzma curl git r-base \
    libboost-dev libcurl4-openssl-dev libboost-all-dev libbz2-dev liblzma-dev \
    libpcre2-dev

USER $MAMBA_USER

# Install Shiba, Tosa, and the beta-binomial regression backend with micromamba
RUN micromamba install -y -n base -c conda-forge -c bioconda -c defaults \
    shiba tosa=1.0.0 r-glmmtmb r-jsonlite

# Install styleframe
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --upgrade pip && \
    pip install styleframe==4.2

# Install shiba2sashimi
RUN pip install shiba2sashimi

# Install the checked-out Shiba source
USER root
COPY . /opt/Shiba
RUN ln -sf /opt/Shiba/shiba.py /opt/conda/bin/shiba.py && \
	ln -sf /opt/Shiba/scshiba.py /opt/conda/bin/scshiba.py
USER $MAMBA_USER

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
