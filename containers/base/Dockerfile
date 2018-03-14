# base starting container for NGS580 analysis pipeline containers
FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
# location for apps
ENV PATH="/opt:${PATH}"
# location for misc scripts
RUN mkdir /opt/bin
ENV PATH="/opt/bin:${PATH}"

# location for misc data
RUN mkdir /data

# mount point for for NYULMC phoenix (Singluarity)
RUN mkdir /ifs
