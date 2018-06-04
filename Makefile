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
TIMESTAMP:=$(shell date +%s)

.PHONY: containers annovar_db ref

# no default action
none:

# ~~~~~ SETUP PIPELINE ~~~~~ #
# install Nextflow in the current directory
./nextflow:
	if [ "$$( module > /dev/null 2>&1; echo $$?)" -eq 0 ]; then module unload java && module load java/1.8 ; fi ; \
	export NXF_VER="$(NXF_VER)" && \
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

update: ./nextflow
	./nextflow self-update

# setup reference data
ref: install
	if [ ! -d "$(REFDIR)" ] ; then echo ">>> system ref dir doesnt exist, setting up local ref dir..." ; \
	./nextflow run ref.nf -profile ref ; fi

# demo pipeline dataset for testing
NGS580-demo-data:
	git clone https://github.com/NYU-Molecular-Pathology/NGS580-demo-data.git

# set up for demo with 'small' dataset (more data output)
demo-small: NGS580-demo-data
	./generate-samplesheets.py NGS580-demo-data/small/ && \
	mv targets.bed "targets.bed.$(TIMESTAMP)" && \
	/bin/cp NGS580-demo-data/small/targets.bed .
	./update-samplesheets.py --tumor-normal-sheet NGS580-demo-data/small/samples.tumor.normal.csv

# set up for demo with 'tiny' dataset (faster processing, ~1hr on NYU HPC)
demo: NGS580-demo-data
	./generate-samplesheets.py NGS580-demo-data/tiny/fastq/ && \
	mv targets.bed "targets.bed.$(TIMESTAMP)" && \
	/bin/cp NGS580-demo-data/tiny/targets.bed .
	./update-samplesheets.py --tumor-normal-sheet NGS580-demo-data/tiny/samples.pairs.csv \
	--pairs-tumor-colname '#SAMPLE-T' \
	--pairs-normal-colname '#SAMPLE-N'

# setup ANNOVAR reference databases
annovar_db: install
	if [ ! -d "$(ANNOVAR_DB_DIR)" ] ; then echo ">>> system ANNOVAR db dir does not exist, setting up local dir..." ;  ./nextflow run annovar_db.nf -profile annovar_db ; fi

annovar_db_power: install
	if [ ! -d "$(ANNOVAR_DB_DIR)" ] ; then echo ">>> system ANNOVAR db dir does not exist, setting up local dir..." ;  ./nextflow run annovar_db.nf -profile annovar_db_conda ; fi

# main setup commands to use
setup: install ref annovar_db
setup-power: install ref annovar_db_power

# prepares a new directory with a copy of this repo to start analysis
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
	@module unload java && module load java/1.8 && \
	log_file="nextflow.$(TIMESTAMP).stdout.log" ; \
	fail_file="failed.$(TIMESTAMP).tsv" ; \
	echo ">>> Running Nextflow with stdout log file: $${log_file}" ; \
	./nextflow run main.nf -profile phoenix -resume -with-dag flowchart-NGS580.dot $(EP) | tee -a "$${log_file}" ; \
	grep '^FAIL' "$${log_file}" > "$${fail_file}" ; \
	echo ">>> Nextflow completed, stdout log file: $${log_file}, failed file: $${fail_file}"


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
	echo 'make run-phoenix-qsub EP="$(EP)"' | qsub -wd "$$PWD" -o :$${qsub_logdir}/ -e :$${qsub_logdir}/ -j y -N "$$job_name" -q all.q 

# parent Nextflow process to be run as a qsub job on phoenix
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

# issue an interupt signal to a process running on a remote server
# e.g. Nextflow running in a qsub job on a compute node
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

clean-old-stdout-logs:
	find . -maxdepth 1 -mindepth 1 -type f -name "nextflow.*.stdout.log" | sort -r | tail -n +2 | xargs rm -f

clean-old-failed-logs:
	find . -maxdepth 1 -mindepth 1 -type f -name "failed.*.tsv" | sort -r | tail -n +2 | xargs rm -f

# clean all files produced by previous pipeline runs
clean: clean-logs clean-traces clean-reports clean-flowcharts clean-old-stdout-logs clean-old-failed-logs

# clean all files produced by all pipeline runs
clean-all: clean clean-output clean-work
	[ -d .nextflow ] && mv .nextflow .nextflowold && rm -rf .nextflowold &
	rm -f .nextflow.log
	rm -f *.png
	rm -f trace*.txt*
	rm -f *.html*
	rm -f flowchart*.dot
	rm -f nextflow.*.stdout.log



# ~~~~~ FINALIZE ~~~~~ #
# steps for finalizing the Nextflow pipeline output 'work' directory
#  Work dir contains subdirectories with intermediary files created when executing pipeline
#  along with command execution logs, etc
#  logs should be saved, along with a record of the contents of each subdir
#  files in the 'publishDir' are often symlinked to items in the 'work' dir;
#  these symlinks should be replaced with copies of the original file
#  all files in the 'work' dir should be removed, except for the log files
# 
#  Makefile configured for parallel processing of files
# 
#  run with `make <recipe> -j8`
# 

LSFILEREGEX:=\.ls\.txt
# files in work subdirs to keep
NXFWORKFILES:='.command.begin|.command.err|.command.log|.command.out|.command.run|.command.sh|.command.stub|.command.trace|.exitcode|$(LSFILE)'
NXFWORKFILESREGEX:='.*\.command\.begin\|.*\.command\.err\|.*\.command\.log\|.*\.command\.out\|.*\.command\.run\|.*\.command\.sh\|.*\.command\.stub\|.*\.command\.trace\|.*\.exitcode\|.*$(LSFILEREGEX)'
# Nextflow 'trace' file with record of completed pipeline tasks

# 

# .PHONY: $(publishDirLinks) $(NXFWORKSUBDIRS)


# Nextflow "publishDir" directory of files to keep
publishDir:=output
# symlinks in the publishDir to convert to files
publishDirLinks:=
FIND_publishDirLinks:=
ifneq ($(FIND_FILES),)
publishDirLinks:=$(shell find $(publishDir)/ -type l)
endif
finalize-output:
	$(MAKE) finalize-output-recurse FIND_publishDirLinks=1
finalize-output-recurse: $(publishDirLinks)
# convert all symlinks to their linked items
$(publishDirLinks):
	@ { \
	destination="$@"; \
	sourcepath="$$(python -c 'import os; print(os.path.realpath("$@"))')" ; \
	if [ ! -e "$${sourcepath}" ]; then echo "ERROR: Source does not exist: $${sourcepath}"; \
	elif [ -f "$${sourcepath}" ]; then rsync -va "$$sourcepath" "$$destination" ; \
	elif [ -d "$${sourcepath}" ]; then { \
	timestamp="$$(date +%s)" ; \
	tmpdir="$${destination}.$${timestamp}" ; \
	rsync -va "$${sourcepath}/" "$${tmpdir}" && \
	rm -f "$${destination}" && \
	mv "$${tmpdir}" "$${destination}" ; } ; \
	fi ; }
.PHONY: $(publishDirLinks)


# Nextflow "work" directory of items to be removed
workDir:=work
# subdirs in the 'work' dir
NXFWORKSUBDIRS:=
FIND_NXFWORKSUBDIRS:=
# regex from the hashes of tasks in the tracefile to match against work subdirs
HASHPATTERN:=
TRACEFILE:=trace-NGS580.txt
ifneq ($(FIND_NXFWORKSUBDIRS),)
NXFWORKSUBDIRS:=$(shell find "$(workDir)/" -maxdepth 2 -mindepth 2) 
HASHPATTERN:=$(shell python -c 'import csv; reader = csv.DictReader(open("$(TRACEFILE)"), delimiter = "\t"); print("|".join([row["hash"] for row in reader]))')
endif
# file to write 'ls' contents of 'work' subdirs to
LSFILE:=.ls.txt
finalize-work-ls:
	$(MAKE) finalize-work-ls-recurse FIND_NXFWORKSUBDIRS=1
finalize-work-ls-recurse: $(NXFWORKSUBDIRS)
# print the 'ls' contents of each subdir to a file, or delete the subdir 
$(NXFWORKSUBDIRS): 
	@if [ "$$(echo '$@' | grep -q -E "$(HASHPATTERN)"; echo $$? )" -eq 0 ]; then \
	ls_file="$@/$(LSFILE)" ; \
	ls -1 "$@" > "$${ls_file}" ; \
	else \
	echo "$@ is not in the trace file, deleting..." ; rm -rf "$@" ; \
	fi
.PHONY: $(NXFWORKSUBDIRS)



# replace all files in 'work' dirs with empty file stubs
NXFWORKFILES:=
FIND_NXFWORKFILES:=
# files in work subdirs to keep
LSFILEREGEX:=\.ls\.txt
NXFWORKFILES:='.command.begin|.command.err|.command.log|.command.out|.command.run|.command.sh|.command.stub|.command.trace|.exitcode|$(LSFILE)'
NXFWORKFILESREGEX:='.*\.command\.begin\|.*\.command\.err\|.*\.command\.log\|.*\.command\.out\|.*\.command\.run\|.*\.command\.sh\|.*\.command\.stub\|.*\.command\.trace\|.*\.exitcode\|.*$(LSFILEREGEX)'
ifneq ($(FIND_NXFWORKFILES),)
NXFWORKFILES:=$(shell find -P "$(workDir)/" -type f ! -regex $(NXFWORKFILESREGEX))
endif
finalize-work-rm:
	$(MAKE) finalize-work-rm-recurse FIND_NXFWORKFILES=1
finalize-work-rm-recurse: $(NXFWORKFILES)
$(NXFWORKFILES):
	@printf 'Creating file stub: $@\n' && rm -f "$@" && touch "$@"
.PHONY: $(NXFWORKFILES)


# remove all symlinks in 'work' dirs
NXFWORKLINKS:=
FIND_NXFWORKLINKS:=
ifneq ($(FIND_NXFWORKLINKS),)
NXFWORKLINKS:=$(shell find "$(workDir)/" -type l)
endif
finalize-work-unlink:
	$(MAKE) finalize-work-unlink-recurse FIND_NXFWORKLINKS=1
finalize-work-unlink-recurse: $(NXFWORKLINKS)
$(NXFWORKLINKS): 
	@printf "Removing symlink: $@\n" && unlink "$@"
.PHONY: $(NXFWORKLINKS)