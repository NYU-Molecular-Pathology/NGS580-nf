FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
unzip \
default-jdk \
xvfb \
xorg \
python

RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.10.zip && \
unzip IGV_2.4.10.zip && \
rm -f unzip IGV_2.4.10.zip

ENV PATH="/IGV_2.4.10/:${PATH}"

# make a dummy batch script in order to load the hg19 genome into the container
# https://software.broadinstitute.org/software/igv/PortCommands
RUN printf 'new\ngenome hg19\nexit\n' > /genome.bat
RUN xvfb-run --auto-servernum --server-num=1 igv.sh -b /genome.bat
