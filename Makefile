# Makefile to run the pipeline
SHELL:=/bin/bash
REFDIR:=/ifs/data/sequence/results/external/NYU/snuderllab/ref
ANNOVAR_DB_DIR:=/ifs/data/molecpathlab/bin/annovar/db/hg19
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

annovar_db_power: install
	if [ ! -d "$(ANNOVAR_DB_DIR)" ] ; then echo ">>> system ANNOVAR db dir does not exist, setting up local dir..." ;  ./nextflow run annovar_db.nf -profile annovar_db_conda ; fi

# main setup commands to use
setup: install ref annovar_db
setup-power: install ref annovar_db_power

# set up a new sequencing directory with a copy of this repo for analysis
deploy:
	@[ -z "$(PROJECT)" ] && printf "invalid PROJECT specified: $(PROJECT)\n" && exit 1 || :
	@[ ! -d "$(SEQDIR)/$(PROJECT)" ] && printf "PROJECT is not a valid location: $(SEQDIR)/$(PROJECT)\n" && exit 1 || :
	@[ -z "$(FASTQDIR)" ] && printf "invalid FASTQDIR specified: $(FASTQDIR)\n" && exit 1 || :
	@[ ! -d "$(FASTQDIR)" ] && printf "FASTQDIR is not a valid directory: $(FASTQDIR)\n" && exit 1 || :
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
# run on phoenix in current session
run-phoenix: install
	module unload java && module load java/1.8 && \
	./nextflow run main.nf -profile phoenix -resume -with-dag flowchart-NGS580.dot $(EP)

# run locally default settings
run-local: install
	./nextflow run main.nf -profile local -resume -with-dag flowchart-NGS580.dot $(EP)

# run on Power server
run-power: install
	source /shared/miniconda2/bin/activate /shared/biobuilds-2017.11 && \
	./nextflow run main.nf -profile power -resume -with-dag flowchart-NGS580.dot $(EP)

# compile flow chart
flowchart:
	[ -f flowchart-NGS580.dot ] && dot flowchart-NGS580.dot -Tpng -o flowchart-NGS580.png || echo "file flowchart-NGS580.dot not present"



# submit the parent Nextflow process to phoenix HPC as a qsub job
submit-phoenix:
	@qsub_logdir="logs" ; \
	mkdir -p "$${qsub_logdir}" ; \
	job_name="NGS580-nf" ; \
	echo 'make run-phoenix-qsub EP=$(EP)' | qsub -wd $$PWD -o :$${qsub_logdir}/ -e :$${qsub_logdir}/ -j y -N "$$job_name" -q all.q 

# parent Nextflow process to be run as a qsub job
run-phoenix-qsub: install
	@output_file="pid.txt" ; \
	module unload java && module load java/1.8 && \
	JOBINFO="$${JOB_ID:-none}\t$${JOB_NAME:-none}\t$${HOSTNAME:-none}\t$${USER:-none}" ; \
	./nextflow run main.nf -profile phoenix -resume -with-dag flowchart-NGS580.dot $(EP) & \
	pid="$$!" ; \
	INFOSTR="$${pid}\t$${JOBINFO}\t$$(date +%s)" ; \
	printf "$${INFOSTR}\n" ; \
	printf "$${INFOSTR}\n" >> $${output_file} ; \
	wait $${pid}

# issue an interupt signal to a process (e.g. Nextflow) running on a remote server
REMOTE:=
PID:=
remote-kill:
	@[ -z "$(REMOTE)" ] && printf "invalid REMOTE specified: $(PROJECT)\n" && exit 1 || :
	@[ -z "$(PID)" ] && printf "invalid PID specified: $(PID)\n" && exit 1 || :
	ssh "$(REMOTE)" 'kill $(PID)'



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
