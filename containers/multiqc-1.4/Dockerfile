# NGS580 container for MultiQC
FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
ADD multiqc.requirements.txt /multiqc.requirements.txt
RUN apt-get update && \
apt-get install -y python \
python-pip && \
pip install --upgrade pip

# ~~~~~ MULTIQC ~~~~~ #
RUN pip install -r /multiqc.requirements.txt
