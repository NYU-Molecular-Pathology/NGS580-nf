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

