# NGS580 container for BWA
FROM debian:jessie

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
# location for apps and data
RUN mkdir /opt/bin
ENV PATH="/opt/bin:${PATH}"

RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y bzip2
RUN apt-get install -y gcc
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y make


# ~~~~~ BWA ~~~~~ #
RUN cd /opt/bin && wget https://superb-dca2.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2 && tar xvf bwa-0.7.17.tar.bz2 && rm -f bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make
ENV PATH="/opt/bin/bwa-0.7.17:${PATH}"
