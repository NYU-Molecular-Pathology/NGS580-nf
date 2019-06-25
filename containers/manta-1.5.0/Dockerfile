FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

RUN apt-get update && \
apt-get install -y wget \
python \
bzip2

RUN cd /opt/ && \
wget https://github.com/Illumina/manta/releases/download/v1.5.0/manta-1.5.0.centos6_x86_64.tar.bz2 && \
tar xvjf manta-1.5.0.centos6_x86_64.tar.bz2 && \
rm -f manta-1.5.0.centos6_x86_64.tar.bz2

ENV PATH="/opt/manta-1.5.0.centos6_x86_64/bin:${PATH}"
