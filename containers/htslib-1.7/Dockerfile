FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

RUN apt-get update && \
apt-get install -y wget \
bzip2

RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"

RUN conda config --add channels defaults
RUN conda config --add channels bioconda

RUN conda install -y -c bioconda -c conda-forge htslib=1.7
