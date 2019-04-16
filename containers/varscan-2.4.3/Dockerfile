# NGS580 container for making reports with R 3.4.3
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

# ~~~~~ VARSCAN ~~~~~ #
RUN conda install -y -c bioconda \
varscan=2.4.3 \
samtools=1.9

# wrapper script for pipeline compatibility
RUN echo '#!/bin/bash' > /usr/varscan.sh
RUN echo 'java -jar /conda/share/varscan-2.4.3-2/VarScan.jar "$@"' >> /usr/varscan.sh
RUN chmod +x /usr/varscan.sh
ENV PATH="/usr:${PATH}"
