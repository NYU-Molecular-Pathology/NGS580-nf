# NGS580 container for R 3.3.2
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
bzip2 \
make \
gcc \
g++ \
gfortran \
libxml2-dev \
libcairo2-dev

# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"

# ~~~~~ R ~~~~~ #
RUN conda install -y -c r r-base \
r=3.3.2 \
r-base=3.3.2 \
pandoc==1.15.0.6 \
r-reshape2 \
r-knitr \
r-rmarkdown \
r-ggplot2

ADD install.R /install.R
RUN chmod +x /install.R
RUN Rscript install.R
