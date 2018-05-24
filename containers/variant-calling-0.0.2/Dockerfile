# NGS580 container for variant calling and downstream processing
FROM broadinstitute/gatk3:3.8-0

# ~~~~~ BASIC SETUP ~~~~~ #
# mount point for for NYULMC phoenix (Singluarity)
RUN mkdir /ifs

# location for misc data
RUN mkdir /data

# location for misc scripts
RUN mkdir /opt/bin
ENV PATH="/opt/bin:/opt:${PATH}"

# GATK shortcut script
ADD gatk.sh /usr/gatk.sh
RUN chmod +x /usr/gatk.sh
ENV PATH="/usr:${PATH}"

RUN apt-get update && \
\
apt-get install -y wget \
bzip2 \
libcairo2-dev

# ~~~~~ MINICONDA ~~~~~ #
# Miniconda2 required for LoFreq 2.1.2
RUN wget https://repo.continuum.io/miniconda/Miniconda2-4.4.10-Linux-x86_64.sh && \
bash Miniconda2-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda2-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# ~~~~~ SAMTOOLS, etc. ~~~~~ #
RUN conda install -y -c bioconda -c r \
samtools==1.3 \
bcftools==1.3 \
htslib=1.7 \
vcflib=1.0.0 \
lofreq=2.1.2 \
r=3.3.2 \
r-base=3.3.2 \
r-ggplot2==2.2.1
# comes with R 3.1.1 by default

# ~~~~~ ANNOVAR ~~~~~ #
RUN cd /opt/ && \
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision150617.tar.gz && \
tar -vzxf annovar.revision150617.tar.gz && \
rm -f annovar.revision150617.tar.gz
# location to mount reference database dir:
RUN mkdir /opt/annovar/db
# easy script for annotating:
ENV PATH="/opt/annovar:${PATH}"

# ~~~~~ R Libs ~~~~~ #
RUN Rscript -e "install.packages(c('curl', 'ggplot2', 'gsalib', 'gplots', 'reshape', 'plyr'), repos='http://cran.us.r-project.org', dependencies = TRUE)"
