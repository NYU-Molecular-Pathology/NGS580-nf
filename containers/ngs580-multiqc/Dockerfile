# NGS580 container for MultiQC
FROM stevekm/ngs580-base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
ADD multiqc.requirements.txt /multiqc.requirements.txt
RUN apt-get update
RUN apt-get install -y python
RUN apt-get install -y python-pip
RUN pip install --upgrade pip

# ~~~~~ MULTIQC ~~~~~ #
RUN pip install -r /multiqc.requirements.txt
