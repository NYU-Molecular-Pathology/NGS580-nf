.PHONY: annovar-150617
none:

test:
	source /shared/miniconda2/bin/activate /shared/biobuilds-2017.11 && \
	conda --version

build:
	source /shared/miniconda2/bin/activate /shared/biobuilds-2017.11 && \
	conda-build


list-env:
	source /shared/miniconda2/bin/activate && \
	conda env list

list:
	source /shared/miniconda2/bin/activate && \
	conda list
	

# create an env to use for building packages
create-conda-build:
	/shared/miniconda2/bin/conda create -n conda-build-env conda-build

# create ANNOVAR env from local package
annovar-150617:
	source /shared/miniconda2/bin/activate conda-build-env && \
	conda-build annovar-150617 && \
	conda create --name annovar-150617 && \
	conda install --use-local annovar-150617 -n annovar-150617
	

list-conda-build:
	source /shared/miniconda2/bin/activate conda-build-env && \
	conda list


# 	conda create -c system -n annovar-150617-env annovar-150617
# 	conda install --use-local annovar-150617 -n annovar-150617
# wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision150617.tar.gz && \
# 	tar xvfz annovar.revision150617.tar.gz && \
# 	rm -f annovar.revision150617.tar.gz
