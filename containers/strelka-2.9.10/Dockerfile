FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

RUN apt-get update && \
apt-get install -y wget \
python \
bzip2

RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 && \
tar xvjf strelka-2.9.10.centos6_x86_64.tar.bz2 && \
rm -f strelka-2.9.10.centos6_x86_64.tar.bz2
