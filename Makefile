# Makefile to run the pipeline
SHELL:=/bin/bash
REFDIR:=/ifs/data/sequence/results/external/NYU/snuderllab/ref
ANNOVAR_DB_DIR:=/ifs/data/molecpathlab/bin/annovar/db/hg19
ANNOVAR_PROTOCOL:=$(shell head -1 annovar_protocol.txt)
ANNOVAR_BUILD_VERSION:=hg19
NXF_VER:=0.29.0
# extra params to pass for Nextflow in some recipes
EP:=
# sequencing run name for deployment
PROJECT:=
# location of production sequencing directory
SEQDIR:=/ifs/data/molecpathlab/quicksilver
# location of production deployment for analysis
PRODDIR:=/ifs/data/molecpathlab/production/NGS580
# location of production demultiplexing for deployment
FASTQDIR:=

.PHONY: containers annovar_db ref

# no default action
none:

# ~~~~~ SETUP PIPELINE ~~~~~ #
# install Nextflow in the current directory
./nextflow:
	export NXF_VER="$(NXF_VER)" && \
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

update: ./nextflow
	./nextflow self-update

# download or symlink ref dir
ref: install
	if [ ! -d "$(REFDIR)" ] ; then echo ">>> system ref dir doesnt exist, setting up local ref dir..." ; \
	./nextflow run ref.nf -profile ref ; fi

# demo pipeline dataset for testing
NGS580-demo-data:
	git clone https://github.com/NYU-Molecular-Pathology/NGS580-demo-data.git

samples.analysis.tsv: NGS580-demo-data
	./generate-samplesheets.py NGS580-demo-data/tiny/fastq/ && \
	mv targets.bed targets.bed.old && \
	/bin/cp NGS580-demo-data/tiny/targets.bed .

demo: samples.analysis.tsv

annovar_db: install
	if [ ! -d "$(ANNOVAR_DB_DIR)" ] ; then echo ">>> system ANNOVAR db dir does not exist, setting up local dir..." ;  ./nextflow run annovar_db.nf -profile annovar_db ; fi

# main setup commands to use
setup: install ref annovar_db

# setup commands needed for NYU phoenix HPC
setup-phoenix: install ref annovar_db


# set up a new sequencing directory with a copy of this repo for analysis
deploy:
	[ -z "$(PROJECT)" ] && printf "invalid PROJECT specified: $(PROJECT)\n" && exit 1 || :
	[ ! -d "$(SEQDIR)/$(PROJECT)" ] && printf "PROJECT is not a valid location: $(SEQDIR)/$(PROJECT)\n" && exit 1 || :
	[ -z "$(FASTQDIR)" ] && printf "invalid FASTQDIR specified: $(FASTQDIR)\n" && exit 1 || :
	[ ! -d "$(FASTQDIR)" ] && printf "FASTQDIR is not a valid directory: $(FASTQDIR)\n" && exit 1 || :
	repo_dir="$${PWD}" && \
	output_dir="$(PRODDIR)/$(PROJECT)/$$(date +"%Y-%m-%d_%H-%M-%S")" && \
	echo ">>> Setting up repo in location: $${output_dir}" && \
	git clone --recursive "$${repo_dir}" "$${output_dir}" && \
	cd "$${output_dir}" && \
	echo ">>> Linking input directory: $(FASTQDIR)" && \
	( mkdir input ; cd input ; ln -s "$(FASTQDIR)" ) && \
	echo ">>> Creating input fastq sheet" && \
	python generate-samplesheets.py "$(FASTQDIR)" && \
	run_cmd="make run-phoenix" && \
	printf ">>> please run the following command to start analysis:\n\n%s\n%s\n" "cd $${output_dir}" "$${run_cmd}"



# ~~~~~ RUN PIPELINE ~~~~~ #
# run on phoenix default settings
run-phoenix: install
	module unload java && module load java/1.8 && \
	./nextflow run main.nf -profile standard -resume -with-dag flowchart-NGS580.dot $(EP)

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
