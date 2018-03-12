# NGS580 container for FastQC
FROM debian:jessie

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
# location for apps and data
RUN mkdir /opt/bin
ENV PATH="/opt/bin:${PATH}"

RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y unzip
RUN apt-get install -y perl-modules

# ~~~~~ JAVA ~~~~~ #
RUN \
    echo "===> add webupd8 repository..."  && \
    echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" | tee /etc/apt/sources.list.d/webupd8team-java.list  && \
    echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" | tee -a /etc/apt/sources.list.d/webupd8team-java.list  && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys EEA14886  && \
    apt-get update  && \
    \
    \
    echo "===> install Java"  && \
    echo debconf shared/accepted-oracle-license-v1-1 select true | debconf-set-selections  && \
    echo debconf shared/accepted-oracle-license-v1-1 seen true | debconf-set-selections  && \
    DEBIAN_FRONTEND=noninteractive  apt-get install -y --force-yes oracle-java8-installer oracle-java8-set-default  && \
    \
    \
    echo "===> clean up..."  && \
    rm -rf /var/cache/oracle-jdk8-installer  && \
    apt-get clean  && \
    rm -rf /var/lib/apt/lists/*

# ~~~~~ FastQC ~~~~~ #
RUN cd /opt/bin && wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip && unzip fastqc_v0.11.7.zip && rm -f fastqc_v0.11.7.zip
RUN cd /opt/bin/FastQC && chmod 755 fastqc
ENV PATH="/opt/bin/FastQC:${PATH}"
