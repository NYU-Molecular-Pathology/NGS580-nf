FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
make \
g++ \
zlib1g-dev \
python

# ~~~~~ BEDTOOLS ~~~~~ #
RUN cd /opt && \
wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && \
tar -zxvf bedtools-2.27.1.tar.gz && \
rm -f bedtools-2.27.1.tar.gz && \
cd /opt/bedtools2 && \
make
ENV PATH="/opt/bedtools2/bin:${PATH}"
