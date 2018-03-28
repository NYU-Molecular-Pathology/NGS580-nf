SHELL:=/bin/bash
CONDA:=/shared/miniconda3
.PHONY: annovar-150617
none:

list-env:
	source "$(CONDA)/bin/activate" && \
	conda env list

list:
	source "$(CONDA)/bin/activate" && \
	conda list
	

# ~~~~~~~ CUSTOM ~~~~~~~~~ #
# create ANNOVAR env from local package
annovar:
	source "$(CONDA)/bin/activate" && \
	conda-build annovar-150617 && \
	conda create -y -c local -n annovar-150617 annovar==150617

test-annovar:
	source "$(CONDA)/bin/activate" annovar-150617 && \
	conda list && \
	echo $${PATH} && \
	annotate_variation.pl --version


remove-annovar:
	source "$(CONDA)/bin/activate" && \
	conda clean -ay && \
	conda build purge && \
	mv /home/kellys04/conda-bld /home/kellys04/conda-bld$$(date +%s) && \
	conda remove -y --name annovar-150617 --all && \
	conda remove -y --name annovar --all



# ~~~~~~ FROM BIOBUILDS ~~~~~~ #
bwa: 
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n bwa-0.7.17 bwa==0.7.17 Samtools==1.6.0

bwa-test: 
	source "$(CONDA)/bin/activate" bwa-0.7.17 && \
	which bwa && \
	which samtools 

trimmomatic: 
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n trimmomatic-0.36 Trimmomatic==0.36 

trimmomatic-test:
	source "$(CONDA)/bin/activate" trimmomatic-0.36 && \
	which trimmomatic

sambamba: 
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n sambamba-0.6.6 Sambamba==0.6.6 Samtools==1.6.0

sambamba-test: 
	source "$(CONDA)/bin/activate" sambamba-0.6.6 && \
	which sambamba 

fastqc: 
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n fastqc-0.11.5 FASTQC==0.11.5

fastqc-test: 
	source "$(CONDA)/bin/activate" fastqc-0.11.5 && \
	which fastqc 

delly2: 
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n delly2-0.7.6 delly==0.7.6

delly2-test: 
	source "$(CONDA)/bin/activate" delly2-0.7.6 && \
	which delly

bedtools: 
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n bedtools-2.26.0 bedtools==2.26.0

bedtools-test: 
	source "$(CONDA)/bin/activate" bedtools-2.26.0 && \
	which bedtools

lofreq:
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n lofreq-2.1.2 lofreq_star==2.1.2

lofreq-test: 
	source "$(CONDA)/bin/activate" lofreq-2.1.2 && \
	which lofreq

vcflib: 
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n vcflib-2.1.2 vcflib==1.0.0_rc1_16.05.18

vcflib-test: 
	source "$(CONDA)/bin/activate" vcflib-2.1.2 && \
	which vcf2tsv


R: 
	/bin/cp ../bin/install.R .
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -n r-3.4.2 \
	r==3.4.2 \
	r-ggplot2==2.2.1 \
	Bioconductor==3.6
	source "$(CONDA)/bin/activate" r-3.4.2 && \
	Rscript ./install.R

R-test:
	source "$(CONDA)/bin/activate" r-3.4.2 && \
	R --version && \
	which R

variant-calling-0.0.2: 
	source "$(CONDA)/bin/activate" && \
	conda create -c biobuilds -c local -n variant-calling-0.0.2 \
	Samtools==1.6.0 \
	lofreq_star==2.1.2 \
	annovar==150617 \
	vcflib==1.0.0_rc1_16.05.18 \
	r==3.4.2 \
	r-ggplot2==2.2.1
	source "$(CONDA)/bin/activate" variant-calling-0.0.2 && \
	Rscript -e "install.packages(c('curl', 'gsalib', 'gplots', 'reshape', 'plyr'), repos='http://cran.us.r-project.org', dependencies = TRUE)"


variant-calling-0.0.2-test:
	source "$(CONDA)/bin/activate" variant-calling-0.0.2 && \

# R: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n r-3.3.2 r==3.3.2
# R-test:
# 	source "$(CONDA)/bin/activate" r-3.3.2 && \
# 	R --version && \
# 	which R
# /home/kellys04/.conda/envs/r-3.3.2/lib/R/bin/exec/R: error while loading shared libraries: libiconv.so.2: cannot open shared object file: No such file or directory
