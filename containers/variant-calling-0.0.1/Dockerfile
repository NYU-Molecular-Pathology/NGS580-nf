# NGS580 container for variant calling and downstream processing
FROM broadinstitute/gatk:4.0.1.2
# - samtools Version: 0.1.19-96b5f2294a

# ~~~~~ BASIC SETUP ~~~~~ #
# mount point for for NYULMC phoenix (Singluarity)
RUN mkdir /ifs

# location for misc data
RUN mkdir /data

# location for misc scripts
RUN mkdir /opt/bin
ENV PATH="/opt/bin:${PATH}"

# ~~~~~ LOFREQ ~~~~~ #
RUN cd /opt/ && \
wget https://raw.githubusercontent.com/CSB5/lofreq/master/dist/lofreq_star-2.1.2_linux-x86-64.tgz && \
tar -vzxf lofreq_star-2.1.2_linux-x86-64.tgz && \
rm -f lofreq_star-2.1.2_linux-x86-64.tgz
ENV PATH="/opt/lofreq_star-2.1.2/bin:${PATH}"

# ~~~~~ HTSLIB ~~~~~ #
RUN cd /opt/ && \
wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
tar xvfj htslib-1.7.tar.bz2 && \
rm -f htslib-1.7.tar.bz2 && \
cd /opt/htslib-1.7 && \
./configure && \
make && \
make install
ENV PATH="/opt/htslib-1.7:${PATH}"

# ~~~~~ ANNOVAR ~~~~~ #
RUN cd /opt/ && \
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision150617.tar.gz && \
tar -vzxf annovar.revision150617.tar.gz && \
rm -f annovar.revision150617.tar.gz
# location to mount reference database dir:
RUN mkdir /opt/annovar/db
# easy script for annotating:
ENV PATH="/opt/annovar:${PATH}"

# ~~~~~ vcflib ~~~~~ #
RUN cd /opt/ && \
git clone --recursive https://github.com/ekg/vcflib.git && \
cd vcflib && \
make
ENV PATH="/opt/vcflib/bin:${PATH}"
