# NGS580 container for trimmomatic-0.36
FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
bzip2

# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"
# ~~~~~ CONDA ~~~~~ #
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# ~~~~~ TRIMMOMATIC ~~~~~ #
RUN conda install -y -c bioconda trimmomatic trimmomatic=0.36
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ADD trimmomatic.conda.sh /opt/bin/trimmomatic.sh
