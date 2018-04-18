# NGS580 container for FastQC
FROM stevekm/ngs580-nf:Miniconda3-4.4.10-Linux-x86_64

MAINTAINER Stephen M. Kelly

# ~~~~~ CONDA ~~~~~ #
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# ~~~~~ FASTQC ~~~~~ #
RUN conda install -y -c bioconda fastqc fastqc=0.11.7
