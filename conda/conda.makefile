SHELL:=/bin/bash
CONDA:=/shared/miniconda3
VAR:=
.PHONY: annovar-150617 r-3.4.2 variant-calling-0.0.2 $(VAR)

none:

# ~~~~~ BUILD ~~~~~ #
check: 
	@if [ ! -d "$(VAR)" ]; then echo "ERROR: VAR is not a valid directory; VAR=$(VAR)"; exit 1; fi

# standard env build methods
create: check
	source "$(CONDA)/bin/activate" && \
	conda env create --file "$(VAR)/env.yml" --name "$(VAR)"

remove:
	source "$(CONDA)/bin/activate" && \
	conda env remove -n "$(VAR)"


# custom package builds
annovar-150617:
	source "$(CONDA)/bin/activate" && \
	conda-build annovar-150617 && \
	conda create -y -c local -n annovar-150617 annovar==150617

# custom env builds
r-3.4.2/install.R:
	/bin/cp ../bin/install.R r-3.4.2/

r-3.4.2: r-3.4.2/install.R
	source "$(CONDA)/bin/activate" && \
	conda env create --file r-3.4.2/_env.yml --name r-3.4.2 && \
	conda activate r-3.4.2 && \
	Rscript r-3.4.2/install.R

variant-calling-0.0.2: annovar-150617
	source "$(CONDA)/bin/activate" && \
	conda env create --file variant-calling-0.0.2/_env.yml --name variant-calling-0.0.2 && \
	conda activate variant-calling-0.0.2 && \
	Rscript -e "install.packages(c('curl', 'gsalib', 'gplots', 'reshape', 'plyr'), repos='http://cran.us.r-project.org', dependencies = TRUE)"




# ~~~~~~ INSTALL ~~~~~~ #
MINICONDA_sh:=Miniconda3-4.5.1-Linux-ppc64le.sh
MINICONDA_sh_url:=https://repo.continuum.io/miniconda/$(MINICONDA_sh)
MINICONDA_sh_md5:=454e3b786937eeaa50fb7bee991ac19e
CONDA_INSTALL_DIR:=$(shell echo '$$HOME')/conda
CONDA_ACTIVATE:=$(CONDA_INSTALL_DIR)/bin/activate

$(MINICONDA_sh):
	wget "$(MINICONDA_sh_url)"


dl: $(MINICONDA_sh)

$(CONDA_INSTALL_DIR): dl
	@if [ ! -d "$(CONDA_INSTALL_DIR)" ]; then \
	bash "$(MINICONDA_sh)" -b -p "$(CONDA_INSTALL_DIR)"; \
	else \
	printf "Install dir already exists: %s\nExiting..." "$(CONDA_INSTALL_DIR)"; \
	exit 1; fi


# install conda in the current directory and install the conda-build package to it
install: $(CONDA_INSTALL_DIR)
	source "$(CONDA_ACTIVATE)" && \
	conda install -y conda-build

# deprecated...
# list-env:
# 	source "$(CONDA)/bin/activate" && \
# 	conda env list

# list:
# 	source "$(CONDA)/bin/activate" && \
# 	conda list
	

# ~~~~~~~ CUSTOM ~~~~~~~~~ #
# create ANNOVAR env from local package
# annovar:
# 	source "$(CONDA)/bin/activate" && \
# 	conda-build annovar-150617 && \
# 	conda create -y -c local -n annovar-150617 annovar==150617

# test-annovar:
# 	source "$(CONDA)/bin/activate" annovar-150617 && \
# 	conda list && \
# 	echo $${PATH} && \
# 	annotate_variation.pl --version


# remove-annovar:
# 	source "$(CONDA)/bin/activate" && \
# 	conda clean -ay && \
# 	conda build purge && \
# 	mv /home/kellys04/conda-bld /home/kellys04/conda-bld$$(date +%s) && \
# 	conda remove -y --name annovar-150617 --all && \
# 	conda remove -y --name annovar --all

# # ~~~~~~ FROM BIOBUILDS ~~~~~~ #
# bwa: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n bwa-0.7.17 bwa==0.7.17 Samtools==1.6.0

# bwa-test: 
# 	source "$(CONDA)/bin/activate" bwa-0.7.17 && \
# 	which bwa && \
# 	which samtools 

# trimmomatic: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n trimmomatic-0.36 Trimmomatic==0.36 

# trimmomatic-test:
# 	source "$(CONDA)/bin/activate" trimmomatic-0.36 && \
# 	which trimmomatic

# sambamba: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n sambamba-0.6.6 Sambamba==0.6.6 Samtools==1.6.0

# sambamba-test: 
# 	source "$(CONDA)/bin/activate" sambamba-0.6.6 && \
# 	which sambamba 

# fastqc: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n fastqc-0.11.5 FASTQC==0.11.5

# fastqc-test: 
# 	source "$(CONDA)/bin/activate" fastqc-0.11.5 && \
# 	which fastqc 

# delly2: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n delly2-0.7.6 delly==0.7.6

# delly2-test: 
# 	source "$(CONDA)/bin/activate" delly2-0.7.6 && \
# 	which delly

# bedtools: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n bedtools-2.26.0 bedtools==2.26.0

# bedtools-test: 
# 	source "$(CONDA)/bin/activate" bedtools-2.26.0 && \
# 	which bedtools

# lofreq:
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n lofreq-2.1.2 lofreq_star==2.1.2

# lofreq-test: 
# 	source "$(CONDA)/bin/activate" lofreq-2.1.2 && \
# 	which lofreq

# vcflib: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n vcflib-2.1.2 vcflib==1.0.0_rc1_16.05.18

# vcflib-test: 
# 	source "$(CONDA)/bin/activate" vcflib-2.1.2 && \
# 	which vcf2tsv


# R: 
# 	/bin/cp ../bin/install.R .
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n r-3.4.2 \
# 	r==3.4.2 \
# 	r-ggplot2==2.2.1 \
# 	Bioconductor==3.6
# 	source "$(CONDA)/bin/activate" r-3.4.2 && \
# 	Rscript ./install.R

# R-test:
# 	source "$(CONDA)/bin/activate" r-3.4.2 && \
# 	R --version && \
# 	which R

# variant-calling-0.0.2: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -c local -n variant-calling-0.0.2 \
# 	Samtools==1.6.0 \
# 	lofreq_star==2.1.2 \
# 	annovar==150617 \
# 	vcflib==1.0.0_rc1_16.05.18 \
# 	r==3.4.2 \
# 	r-ggplot2==2.2.1
# 	source "$(CONDA)/bin/activate" variant-calling-0.0.2 && \
# 	Rscript -e "install.packages(c('curl', 'gsalib', 'gplots', 'reshape', 'plyr'), repos='http://cran.us.r-project.org', dependencies = TRUE)"


# variant-calling-0.0.2-test:
# 	source "$(CONDA)/bin/activate" variant-calling-0.0.2 && \

# R: 
# 	source "$(CONDA)/bin/activate" && \
# 	conda create -c biobuilds -n r-3.3.2 r==3.3.2
# R-test:
# 	source "$(CONDA)/bin/activate" r-3.3.2 && \
# 	R --version && \
# 	which R
# /home/kellys04/.conda/envs/r-3.3.2/lib/R/bin/exec/R: error while loading shared libraries: libiconv.so.2: cannot open shared object file: No such file or directory
