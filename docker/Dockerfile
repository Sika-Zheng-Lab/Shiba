###############################################################
# Dockerfile to build container images for Shiba (v0.5.5)
# Based on python 3.11-buster
################################################################

FROM python:3.11-buster

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

# Make directory to put softwares
RUN mkdir /opt_shiba/

# Install R (version 4.1.3)
RUN wget https://cran.r-project.org/src/base/R-4/R-4.1.3.tar.gz && \
    tar -zxvf R-4.1.3.tar.gz && \
    rm -rf R-4.1.3.tar.gz && \
    cd R-4.1.3 && \
	./configure \
    --prefix=/opt_shiba/R/4.1.3 \
	--with-pcre1 \
    --enable-R-shlib \
    --enable-memory-profiling \
    --with-blas \
    --with-lapack && \
	make && \
	make install

ENV PATH /opt_shiba/R/4.1.3/bin:$PATH

# Install DESeq2
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')" && \
	R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos = NULL, type = "source")' && \
    R -e "BiocManager::install('DESeq2')"

# Install python packages
RUN /usr/local/bin/python -m pip install --upgrade pip && \
    pip install pandas==1.5.3 statsmodels==0.13.5 matplotlib==3.6.3 \
	seaborn==0.12.2 numexpr==2.8.4 styleframe==4.1 plotly==5.13.0 dash==2.8.1 \
	pysam==0.21.0 scanpy==1.9.5 \
	numpy==1.26.4 PyYAML==6.0.2

# Install Subread (version 2.0.3)
RUN cd /opt_shiba && \
	wget https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz && \
	tar -zxvf subread-2.0.3-Linux-x86_64.tar.gz && \
	rm -rf subread-2.0.3-Linux-x86_64.tar.gz

# Install StringTie (version 2.2.1)
RUN cd /opt_shiba && \
	wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.tar.gz && \
	tar -zxvf stringtie-2.2.1.tar.gz && \
	rm -rf stringtie-2.2.1.tar.gz && \
	cd stringtie-2.2.1 && \
	make release

# Install Regtools (version 1.0.0)
RUN cd /opt_shiba && \
	wget https://github.com/griffithlab/regtools/archive/refs/tags/1.0.0.tar.gz && \
    tar -zxvf 1.0.0.tar.gz && \
    rm -rf 1.0.0.tar.gz && \
	cd regtools-1.0.0/ && \
	mkdir build && \
	cd build/ && \
	cmake .. && \
	make

# Install Shiba
RUN cd /opt_shiba && \
	git clone https://github.com/Sika-Zheng-Lab/Shiba.git -b v0.5.5

# Install shiba2sashimi
RUN pip install shiba2sashimi==0.1.2

# Set environment variables
ENV PATH $PATH:/opt_shiba/subread-2.0.3-Linux-x86_64/bin:/opt_shiba/stringtie-2.2.1:/opt_shiba/regtools-1.0.0/build:/opt_shiba/Shiba

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
