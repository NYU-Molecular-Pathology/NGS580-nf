FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
bzip2 

# ~~~~~ MINICONDA ~~~~~ #
# Miniconda2 required for LoFreq 2.1.3
RUN wget https://repo.continuum.io/miniconda/Miniconda2-4.4.10-Linux-x86_64.sh && \
bash Miniconda2-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda2-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN conda install -y -c bioconda \
bcftools==1.3 \
lofreq=2.1.3 
