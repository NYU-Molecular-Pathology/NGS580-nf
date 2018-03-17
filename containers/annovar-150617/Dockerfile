# NGS580 container for ANNOVAR
FROM stevekm/ngs580-nf:base

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget

# ~~~~~ ANNOVAR ~~~~~ #
RUN cd /opt/ && \
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision150617.tar.gz && \
tar -vzxf annovar.revision150617.tar.gz && \
rm -f annovar.revision150617.tar.gz
# location to mount reference database dir:
RUN mkdir /opt/annovar/db
# easy script for annotating:
ADD annotate_vcf.sh /opt/annovar/annotate_vcf.sh
ENV PATH="/opt/annovar:${PATH}"
