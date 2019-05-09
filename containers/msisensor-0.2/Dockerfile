# NGS580 container for MSISensor
FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
\
apt-get install -y wget \
bzip2 \
g++ \
# for Samtools
gcc \
make \
libncurses5-dev \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
git

# ~~~~~ SAMTOOLS ~~~~~ #
RUN cd /opt && \
wget https://newcontinuum.dl.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2 && \
tar xvjf samtools-0.1.19.tar.bz2 && \
rm -f samtools-0.1.19.tar.bz2 && \
cd samtools-0.1.19 && \
make
ENV SAMTOOLS_ROOT="/opt/samtools-0.1.19"
ENV PATH="/opt/samtools-0.1.19:${PATH}"

# ~~~~~ MSISENSOR ~~~~~ #
RUN cd /opt && \
git clone https://github.com/ding-lab/msisensor.git && \
cd msisensor && \
make
ENV PATH="/opt/msisensor:${PATH}"
