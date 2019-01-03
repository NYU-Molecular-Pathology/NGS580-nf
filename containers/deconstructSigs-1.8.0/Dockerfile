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
RUN conda config --add channels defaults && \
conda config --add channels conda-forge && \
conda config --add channels bioconda

# ~~~~~ R ~~~~~ #
RUN conda install -c r -c bioconda r=3.4.1 r-base=3.4.1 \
bioconductor-bsgenome=1.48.0 \
bioconductor-bsgenome.hsapiens.ucsc.hg19=1.4.0 \
bioconductor-genomicranges=1.32.7 \
bioconductor-genomeinfodb=1.16.0 \
r-deconstructsigs=1.8.0 \
r-ggplot2=2.2.1
ENV LC_ALL="C.UTF-8"
ENV LANG="C.UTF-8"
