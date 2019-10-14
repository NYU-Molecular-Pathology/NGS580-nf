# NGS580 container for variant calling and downstream processing
FROM broadinstitute/gatk3:3.8-0

# https://unix.stackexchange.com/questions/508724/failed-to-fetch-jessie-backports-repository
RUN echo "deb [check-valid-until=no] http://cdn-fastly.deb.debian.org/debian jessie main" > /etc/apt/sources.list.d/jessie.list
RUN echo "deb [check-valid-until=no] http://archive.debian.org/debian jessie-backports main" > /etc/apt/sources.list.d/jessie-backports.list
RUN sed -i '/deb http:\/\/deb.debian.org\/debian jessie-updates main/d' /etc/apt/sources.list
RUN apt-get -o Acquire::Check-Valid-Until=false update

RUN apt-get install -y wget \
bzip2 \
libcairo2-dev

# ~~~~~ MINICONDA ~~~~~ #
# Miniconda2 required for LoFreq 2.1.2
RUN wget https://repo.continuum.io/miniconda/Miniconda2-4.4.10-Linux-x86_64.sh && \
bash Miniconda2-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda2-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# ~~~~~ SAMTOOLS, etc. ~~~~~ #
RUN conda install -y -c bioconda -c r \
samtools=1.3 \
bcftools=1.3 \
htslib=1.7 \
vcflib=1.0.0 \
lofreq=2.1.2 \
r=3.3.2 \
r-base=3.3.2 \
r-ggplot2=2.2.1 \
r-curl=2.3 \
r-gsalib=2.1 \
r-gplots=3.0.1 \
r-reshape=0.8.6 \
r-plyr=1.8.4
# comes with R 3.1.1 by default

# wrapper script for pipeline compatibility
RUN echo '#!/bin/bash' > /usr/gatk.sh
RUN echo 'java -jar $(dirname $0)/GenomeAnalysisTK.jar "$@"' >> /usr/gatk.sh
RUN chmod +x /usr/gatk.sh
ENV PATH="/usr:${PATH}"
