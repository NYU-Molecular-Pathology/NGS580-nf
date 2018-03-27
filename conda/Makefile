SHELL:= /bin/bash
CONDA:=/shared/miniconda3
.PHONY: annovar-150617
none:

# test:
# 	source "$(CONDA_ACTIVATE)" /shared/biobuilds-2017.11 && \
# 	conda --version

# build:
# 	source "$(CONDA_ACTIVATE)" /shared/biobuilds-2017.11 && \
# 	conda-build


list-env:
	source "$(CONDA)/bin/activate" && \
	conda env list

list:
	source "$(CONDA)/bin/activate" && \
	conda list
	

# create an env to use for building packages
# create-conda-build:
# 	"$(CONDA)/bin/conda" create -n conda-build-env conda-build

# create ANNOVAR env from local package
annovar:
	source "$(CONDA)/bin/activate" && \
	package_path="$$(conda-build annovar-150617 --output --no-force-upload --no-anaconda-upload)" && \
	echo "$$package_path" && \
	conda install --use-local "$$package_path" -n annovar-150617


test-annovar:
	source "$(CONDA)/bin/activate" annovar-150617 && \
	conda list && \
	echo $${PATH} && \
	echo $${PYTHONPATH}

remove-annovar:
	source "$(CONDA)/bin/activate" && \
	conda remove --name annovar-150617 --all

# conda create --name annovar-150617 && \
# conda install --use-local annovar-150617 -n annovar-150617


# test:
# 	output="$$(ls /foo)" && \
# 	echo "$$output_path"


# list-conda-build:
# 	source "$(CONDA)/bin/activate" conda-build-env && \
# 	conda list


# 	conda create -c system -n annovar-150617-env annovar-150617
# 	conda install --use-local annovar-150617 -n annovar-150617
# wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision150617.tar.gz && \
# 	tar xvfz annovar.revision150617.tar.gz && \
# 	rm -f annovar.revision150617.tar.gz
