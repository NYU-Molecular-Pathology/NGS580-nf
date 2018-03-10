FROM stevekm/ngs580-base

MAINTAINER Stephen M. Kelly

# ~~~~~ BWA ~~~~~ #
RUN cd /opt/bin && wget https://superb-dca2.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2 && tar xvf bwa-0.7.17.tar.bz2 && rm -f bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make
ENV PATH="/opt/bin/bwa-0.7.17:${PATH}"
