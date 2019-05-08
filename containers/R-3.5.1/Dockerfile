FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

RUN apt-get update && \
apt-get install -y wget \
bzip2

RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.3.31-Linux-x86_64.sh && \
bash Miniconda3-4.3.31-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda3-4.3.31-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"

RUN conda install -y \
-c r \
-c conda-forge \
r-base=3.5.1 \
r=3.5.1 \
r-binom=1.1_1

ENV LC_ALL="C.UTF-8"
ENV LANG="C.UTF-8"
