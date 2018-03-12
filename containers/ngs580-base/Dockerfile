# base starting container for NGS580 analysis pipeline containers
FROM debian:jessie

MAINTAINER Stephen M. Kelly

# location for apps and data
RUN mkdir /opt/bin
ENV PATH="/opt/bin:${PATH}"

RUN apt-get update
RUN apt-get install -y wget
# RUN apt-get install -y git
RUN apt-get install -y unzip
RUN apt-get install -y make
RUN apt-get install -y bzip2
RUN apt-get install -y gcc
# RUN apt-get install -y libncurses5-dev
RUN apt-get install -y zlib1g-dev
# RUN apt-get install -y g++
# RUN apt-get install -y libbz2-dev
# RUN apt-get install -y liblzma-dev
