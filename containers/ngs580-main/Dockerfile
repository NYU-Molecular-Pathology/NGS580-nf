FROM stevekm/ngs580-base

MAINTAINER Stephen M. Kelly

RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y git
RUN apt-get install -y unzip
RUN apt-get install -y make
RUN apt-get install -y bzip2
RUN apt-get install -y gcc
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y g++
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev

# ~~~~~ hstlib ~~~~~ #
RUN cd /opt/bin && wget https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && tar xvfj htslib-1.3.1.tar.bz2 && rm -f htslib-1.3.1.tar.bz2
RUN cd /opt/bin/htslib-1.3.1 && ./configure && make

# ~~~~~ SAMTOOLS ~~~~~ #
RUN cd /opt/bin && wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && tar xvfj samtools-1.3.1.tar.bz2 && rm -f samtools-1.3.1.tar.bz2
RUN cd /opt/bin/samtools-1.3.1 && ./configure && make && make install

# ~~~~~ vcflib ~~~~~ #
RUN cd /opt/bin && git clone --recursive https://github.com/ekg/vcflib.git && cd vcflib && make

# ~~~~~ bcftools ~~~~~ #
RUN cd /opt/bin && wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 && tar xvfj bcftools-1.3.1.tar.bz2 && cd bcftools-1.3.1 && make

# ~~~~~ IGV ~~~~~ #
RUN cd /opt/bin && wget http://data.broadinstitute.org/igv/projects/downloads/2.3/IGV_2.3.81.zip -O tmp && mv tmp IGV_2.3.81.zip && unzip IGV_2.3.81.zip && rm -f IGV_2.3.81.zip

# ~~~~~ PANDOC ~~~~~ #
RUN cd /opt/bin && wget https://s3.amazonaws.com/rstudio-buildtools/pandoc-1.13.1.zip && unzip pandoc-1.13.1.zip && rm -f pandoc-1.13.1.zip
