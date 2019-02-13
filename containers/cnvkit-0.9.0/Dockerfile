FROM ubuntu:16.04
MAINTAINER Varshini Vasu
##https://raw.githubusercontent.com/etal/cnvkit/eaa2e5ebfdc09a1cfc495d27cf205d0c37b8059f/docker/Dockerfile

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y r-base-core
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('PSCBS', 'cghFLasso')"

RUN apt-get install -y \
    liblzma-dev \
    python-biopython \
    python-dev \
    python-matplotlib \
    python-numpy \
    python-pip \
    python-reportlab \
    python-scipy \
    python-tk \
    zlib1g-dev

RUN pip install \
future==0.17.1 \
futures==3.2.0 \
pyfaidx==0.5.5.2 \
pysam==0.15.2 \
pandas==0.20.3 \
cnvkit==0.9.0

# Let matplotlib build its font cache
RUN cnvkit.py version
