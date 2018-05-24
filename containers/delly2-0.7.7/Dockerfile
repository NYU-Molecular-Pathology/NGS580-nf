# NGS580 container for Delly2
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
\
apt-get install -y git \
wget \
make \
gcc \
g++ \
libncurses5-dev \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
bzip2

# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda2-4.4.10-Linux-x86_64.sh && \
bash Miniconda2-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda2-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# # ~~~~~ DELLY2, etc. ~~~~~ #
RUN conda install -y -c bioconda \
delly==0.7.7
