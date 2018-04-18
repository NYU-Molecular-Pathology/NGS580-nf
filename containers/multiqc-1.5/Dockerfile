# NGS580 container for multiqc-1.5
FROM stevekm/ngs580-nf:Miniconda3-4.4.10-Linux-x86_64

MAINTAINER Stephen M. Kelly

# ~~~~~ CONDA ~~~~~ #
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# ~~~~~ MULTIQC ~~~~~ #
RUN conda install -y -c bioconda multiqc multiqc=1.5
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
