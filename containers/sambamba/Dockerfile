# NGS580 container for Sambamba and Samtools
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update
# for Sambamaba
RUN apt-get install -y wget
RUN apt-get install -y bzip2
# for Samtools
RUN apt-get install -y gcc
RUN apt-get install -y make
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev


# ~~~~~ SAMBAMBA ~~~~~ #
RUN cd /opt/bin && wget https://github.com/biod/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2 && tar xvjf sambamba_v0.6.6_linux.tar.bz2 && rm -f sambamba_v0.6.6_linux.tar.bz2 && ln -s sambamba_v0.6.6 sambamba

# ~~~~~ SAMTOOLS ~~~~~ #
RUN cd /opt/bin && wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && tar xvjf samtools-1.7.tar.bz2 && rm -f samtools-1.7.tar.bz2
RUN cd /opt/bin/samtools-1.7/ && ./configure && make && make install
