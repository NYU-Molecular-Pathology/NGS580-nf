# NGS580 container for MultiQC
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
ADD multiqc.requirements.txt /multiqc.requirements.txt
RUN apt-get update && \
apt-get install -y python \
python-pip && \
pip install --upgrade pip

# ~~~~~ MULTIQC ~~~~~ #
RUN pip install -r /multiqc.requirements.txt
