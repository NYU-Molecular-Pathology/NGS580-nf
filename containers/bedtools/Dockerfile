# NGS580 container for bedtools
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y make
RUN apt-get install -y g++
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y python


# ~~~~~ BEDTOOLS ~~~~~ #
RUN cd /opt/bin && wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && tar -zxvf bedtools-2.26.0.tar.gz && rm -f bedtools-2.26.0.tar.gz
RUN cd /opt/bin/bedtools2 && make
ENV PATH="/opt/bin/bedtools2/bin:${PATH}"
