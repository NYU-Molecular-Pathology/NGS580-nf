# NGS580 container for Delly2
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
\
apt-get install -y git \
wget \
make \
gcc \
g++ \
libncurses5-dev \
zlib1g-dev \
libbz2-dev \
liblzma-dev



# ~~~~~ DELLY2 ~~~~~ #
RUN cd /opt && \
git clone --recursive https://github.com/dellytools/delly.git && \
cd delly && \
git checkout tags/v0.7.7 && \
git submodule update --recursive && \
make
ENV PATH="/opt/delly/src:${PATH}"
