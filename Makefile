# Makefile to run the pipeline
SHELL:=/bin/bash
REFDIR:=/ifs/data/sequence/results/external/NYU/snuderllab/ref
ANNOVAR_DB_DIR:=/ifs/data/molecpathlab/bin/annovar/db/hg19
ANNOVAR_PROTOCOL:=$(shell head -1 annovar_protocol.txt)
ANNOVAR_BUILD_VERSION:=hg19
NXF_VER:=0.28.0
# extra params to pass for Nextflow in some recipes
EP:=
.PHONY: containers

# no default action
none:

# ~~~~~ SETUP PIPELINE ~~~~~ #
# install Nextflow in the current directory
./nextflow:
	export NXF_VER="$(NXF_VER)" && \
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

# set up MultiQC in a virtualenv
bin/multiqc-venv/bin/activate:
	cd bin && \
	make -f multiqc.makefile setup

multiqc: bin/multiqc-venv/bin/activate

# download or symlink ref dir
ref:
	[ -d "$(REFDIR)" ] && ln -fs $(REFDIR) ref || { wget https://genome.med.nyu.edu/results/external/NYU/snuderllab/ref.tar.gz && \
	tar -vxzf ref.tar.gz && \
	rm -f ref.tar.gz ; }

# clean up ref tarballs 
ref-clean:
	rm -f ref.tar.gz*

# build all Docker containers
build-containers:
	cd containers && make build

# pull all Docker containers from Docker hub
pull-containers:
	cd containers && make pull

# demo pipeline dataset for testing
NGS580-demo-data:
	git clone https://github.com/NYU-Molecular-Pathology/NGS580-demo-data.git

samples.analysis.tsv: NGS580-demo-data
	./generate-samplesheets.py NGS580-demo-data/tiny/fastq/ && \
	mv targets.bed targets.bed.old && \
	/bin/cp NGS580-demo-data/tiny/targets.bed .

demo: samples.analysis.tsv

# set up ANNOVAR reference db dir based on first line in 'annovar_protocol.txt'
annovar_db: bin/annotate_variation.pl
	[ -d "$(ANNOVAR_DB_DIR)" ] && ln -fs $(ANNOVAR_DB_DIR) annovar_db && rm -f annovar.revision*.tar.gz || { \
	mkdir -p annovar_db && \
	for item in $$( echo "$(ANNOVAR_PROTOCOL)" | tr ',' ' ' ) ; do \
	( \
	export PATH="annovar:$${PATH}" ; \
	downdb_param="$$(grep "$$item" annovar_key.tsv | cut -f1)" ; \
	echo "$$downdb_param" ; \
	bin/annotate_variation.pl -downdb -buildver $(ANNOVAR_BUILD_VERSION) -webfrom annovar "$$downdb_param" annovar_db ; \
	) ; \
	done; \
	}

# download the ANNOVAR db's in parrallel
annovar_db_p: bin/annotate_variation.pl
	[ -d "$(ANNOVAR_DB_DIR)" ] && ln -fs $(ANNOVAR_DB_DIR) annovar_db && rm -f annovar.revision*.tar.gz || { \
	mkdir -p annovar_db && \
	for item in $$( echo "$(ANNOVAR_PROTOCOL)" | tr ',' ' ' ) ; do \
	( \
	export PATH="annovar:$${PATH}" ; \
	downdb_param="$$(grep "$$item" annovar_key.tsv | cut -f1)" ; \
	echo "$$downdb_param" ; \
	bin/annotate_variation.pl -downdb -buildver $(ANNOVAR_BUILD_VERSION) -webfrom annovar "$$downdb_param" annovar_db & \
	) ; \
	done; \
	}

# download ANNOVAR; needed to download the reference db's
bin/annotate_variation.pl:
	cd bin && \
	make -f annovar.makefile install

# main setup commands to use
setup: install ref annovar_db 

# setup commands needed for NYU phoenix HPC
setup-phoenix: install ref annovar_db




# ~~~~~ RUN PIPELINE ~~~~~ #
# run on phoenix default settings
run-phoenix: setup-phoenix
	./nextflow run main.nf -profile standard -resume -with-dag flowchart-NGS580.dot $(EP) 

# run on phoenix Singularity head node config
run-phoenix-head: setup-phoenix
	./nextflow run main.nf -profile headnode -with-dag flowchart-NGS580.dot $(EP) 

# run locally default settings
run-local: install ref
	./nextflow run main.nf -profile local -resume -with-dag flowchart-NGS580.dot $(EP) 

# run on Power server
run-power: install ref
	source /shared/miniconda2/bin/activate /shared/biobuilds-2017.11 && \
	./nextflow run main.nf -profile power -resume -with-dag flowchart-NGS580.dot $(EP)

# compile flow chart
flowchart: 
	[ -f flowchart-NGS580.dot ] && dot flowchart-NGS580.dot -Tpng -o flowchart-NGS580.png || echo "file flowchart-NGS580.dot not present"


# ~~~~~ CLEANUP ~~~~~ #
# commands to clean out items in the current directory after running the pipeline
clean-traces:
	rm -f trace*.txt.*

clean-logs:
	rm -f .nextflow.log.*

clean-reports:
	rm -f *.html.*

clean-flowcharts:
	rm -f *.dot.*

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

# clean all files produced by previous pipeline runs
clean: clean-logs clean-traces clean-reports clean-flowcharts

# clean all files produced by all pipeline runs
clean-all: clean clean-output clean-work
	[ -d .nextflow ] && mv .nextflow .nextflowold && rm -rf .nextflowold &
	rm -f .nextflow.log
	rm -f *.png
	rm -f trace*.txt*
	rm -f *.html*
	rm -f flowchart*.dot

# remove all ANNOVAR dirs and files used for ANNOVAR ref db setup
clean-annovar:
	[ -d annovar_db ] && /bin/mv annovar_db annovar_dbold && rm -rf annovar_dbold &
	[ -d annovar ] && /bin/mv annovar annovarold && rm -rf annovarold &
	rm -f annovar.revision*.tar.gz

clean-ref: ref-clean
