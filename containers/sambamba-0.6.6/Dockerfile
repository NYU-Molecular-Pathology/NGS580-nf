# NGS580 container for Sambamba and Samtools
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
\
# for Sambamaba
apt-get install -y wget \
bzip2 \
# for Samtools
gcc \
make \
libncurses5-dev \
zlib1g-dev \
libbz2-dev \
liblzma-dev


# ~~~~~ SAMBAMBA ~~~~~ #
RUN cd /opt && \
wget https://github.com/biod/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2 && \
tar xvjf sambamba_v0.6.6_linux.tar.bz2 && \
rm -f sambamba_v0.6.6_linux.tar.bz2 && \
ln -s sambamba_v0.6.6 sambamba

# ~~~~~ SAMTOOLS ~~~~~ #
RUN cd /opt && \
wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && \
tar xvjf samtools-1.7.tar.bz2 && \
rm -f samtools-1.7.tar.bz2 && \
cd /opt/samtools-1.7/ && \
./configure && \
make && \
make install
