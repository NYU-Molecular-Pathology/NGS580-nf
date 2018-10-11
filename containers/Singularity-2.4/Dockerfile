FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    squashfs-tools \
    wget \
    python \
    libarchive-dev \
    sudo \
    debootstrap

ENV VERSION=2.4
RUN wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz && \
tar xvf singularity-$VERSION.tar.gz && \
cd singularity-$VERSION && \
./configure --prefix=/usr/local && \
make && \
make install
