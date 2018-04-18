# NGS580 container for Trimmomatic
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
unzip

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


# ~~~~~ TRIMMOMATIC ~~~~~ #
RUN cd /opt && \
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
unzip Trimmomatic-0.36.zip && \
rm -f Trimmomatic-0.36.zip && \
chmod +x /opt/Trimmomatic-0.36/trimmomatic-0.36.jar && \
cd /opt/Trimmomatic-0.36 && \
ln -s trimmomatic-0.36.jar trimmomatic.jar
ADD trimmomatic.sh /opt/Trimmomatic-0.36/trimmomatic.sh
RUN chmod +x /opt/Trimmomatic-0.36/trimmomatic.sh
WORKDIR /opt/Trimmomatic-0.36/
ENV PATH="/opt/Trimmomatic-0.36:${PATH}"
