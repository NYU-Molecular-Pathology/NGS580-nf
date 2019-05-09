# NGS580 container for FastQC
FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# https://stackoverflow.com/questions/14066525/java-lang-error-probable-fatal-errorno-fonts-found
RUN apt-get update && \
apt-get install -y ttf-dejavu && \
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

# ~~~~~ FASTQC ~~~~~ #
RUN conda install -y -c bioconda fastqc fastqc=0.11.7
