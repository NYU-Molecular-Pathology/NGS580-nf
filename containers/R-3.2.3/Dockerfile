# NGS580 container for R
FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y --no-install-recommends r-base=3.2.3-4
ADD install.R /install.R
RUN /install.R
