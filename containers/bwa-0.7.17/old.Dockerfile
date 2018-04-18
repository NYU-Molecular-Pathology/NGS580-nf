# NGS580 container for BWA
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
bzip2 \
gcc \
zlib1g-dev \
make


# ~~~~~ BWA ~~~~~ #
RUN cd /opt && \
wget https://superb-dca2.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2 && \
tar xvf bwa-0.7.17.tar.bz2 && \
rm -f bwa-0.7.17.tar.bz2 && \
cd bwa-0.7.17 && \
make
ENV PATH="/opt/bwa-0.7.17:${PATH}"
