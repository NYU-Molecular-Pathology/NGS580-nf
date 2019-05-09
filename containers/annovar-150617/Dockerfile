# NGS580 container for ANNOVAR
FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
perl \
bzip2

# ~~~~~ ANNOVAR ~~~~~ #
RUN cd /opt/ && \
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision150617.tar.gz && \
tar -vzxf annovar.revision150617.tar.gz && \
rm -f annovar.revision150617.tar.gz
# location to mount reference database dir:
RUN mkdir /opt/annovar/db
ENV PATH="/opt/annovar:${PATH}"

# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"

# ~~~~~ R ~~~~~ #
RUN conda install -y -c r r-base r=3.3.2 r-base=3.3.2
