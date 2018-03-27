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
	

# create ANNOVAR env from local package
annovar:
	source "$(CONDA)/bin/activate" && \
	package_path="$$(conda-build annovar-150617 --output --no-force-upload --no-anaconda-upload)" && \
	echo "$$package_path" && \
	conda create -y -c local -n annovar-150617 annovar

# annovar==150617
# conda install --use-local "$$package_path" -n annovar-150617
test-annovar:
	source "$(CONDA)/bin/activate" annovar-150617 && \
	conda list && \
	echo $${PATH} && \
	echo $${PYTHONPATH} && \
	annotate_variation.pl --version


remove-annovar:
	source "$(CONDA)/bin/activate" && \
	conda clean -ay && \
	conda build purge && \
	mv /home/kellys04/conda-bld /home/kellys04/conda-bld$$(date +%s) && \
	conda remove -y --name annovar-150617 --all

