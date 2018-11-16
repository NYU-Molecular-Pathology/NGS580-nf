# Makefile to run the pipeline
SHELL:=/bin/bash
export NXF_VER:=18.10.1
# extra params to pass for Nextflow in some recipes
EP:=
TIMESTAMP:=$(shell date +%s)
REMOTE_ssh:=git@github.com:NYU-Molecular-Pathology/NGS580-nf.git
REMOTE_http:=https://github.com/NYU-Molecular-Pathology/NGS580-nf.git
.PHONY: annovar_db ref

# no default action
none:

# ~~~~~ SETUP PIPELINE ~~~~~ #
HOSTNAME:=$(shell echo $$HOSTNAME)
USER_HOME=$(shell echo "$$HOME")
USER_DATE:=$(shell date +%s)
NXF_FRAMEWORK_DIR:=$(USER_HOME)/.nextflow/framework/$(NXF_VER)
# gets stuck on NFS drive and prevents install command from finishing
remove-framework:
	@if [ -e "$(NXF_FRAMEWORK_DIR)" ]; then \
	new_framework="$(NXF_FRAMEWORK_DIR).$(USER_DATE)" ; \
	echo ">>> Moving old Nextflow framework dir $(NXF_FRAMEWORK_DIR) to $${new_framework}" ; \
	mv "$(NXF_FRAMEWORK_DIR)" "$${new_framework}" ; \
	fi

# install Nextflow in the current directory
# - removes user Nextflow framework dir if exists
# - if on 'phoenix' HPC, need to load Java module
./nextflow:
	@[ -d "$(NXF_FRAMEWORK_DIR)" ] && $(MAKE) remove-framework || :
	@if grep -q 'phoenix' <<<'$(HOSTNAME)'; then module unload java && module load java/1.8; fi ; \
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

nextflow-self-update: ./nextflow
	./nextflow self-update

# setup reference data
# NYUMC Big Purple data directory location
REFDIR:=/gpfs/data/molecpathlab/ref
ref: install
	@if [ ! -d "$(REFDIR)" ] ; then echo ">>> system ref dir doesnt exist, setting up local ref dir..." ; \
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
# NYUMC Big Purple data directory location
ANNOVAR_DB_DIR:=/gpfs/data/molecpathlab/ref/annovar/db
annovar_db: install
	if [ ! -d "$(ANNOVAR_DB_DIR)" ] ; then echo ">>> system ANNOVAR db dir does not exist, setting up local dir..." ;  ./nextflow run annovar_db.nf -profile annovar_db ; fi

annovar_db_power: install
	if [ ! -d "$(ANNOVAR_DB_DIR)" ] ; then echo ">>> system ANNOVAR db dir does not exist, setting up local dir..." ;  ./nextflow run annovar_db.nf -profile annovar_db_conda ; fi

annovar_db_bigpurple: install
	if [ ! -d "$(ANNOVAR_DB_DIR)" ] ; then echo ">>> system ANNOVAR db dir does not exist, setting up local dir..." ;  ./nextflow run annovar_db.nf -profile annovar_db_bigpurple ; fi


# main setup commands to use
setup: install ref annovar_db
setup-power: install ref annovar_db_power

check-runid:
	@[ -z "$(RUNID)" ] && printf "invalid RUNID specified: $(RUNID)\n" && exit 1 || :
check-fastqdir:
	@[ -z "$(FASTQDIR)" ] && printf "invalid FASTQDIR specified: $(FASTQDIR)\n" && exit 1 || :
	@[ ! -d "$(FASTQDIR)" ] && printf "FASTQDIR is not a valid directory: $(FASTQDIR)\n" && exit 1 || :

# prepares a new directory with a copy of this repo to start analysis
# `make deploy RUNID=180316_NB501073_0036_AH3VFKBGX5 FASTQDIR=/ifs/data/molecpathlab/production/Demultiplexing/180316_NB501073_0036_AH3VFKBGX5/output/Unaligned/NS18-7`
# - check that valid args were passed
# - clone current repo into destination for analysis
# - create symlink to fastq dir for input
# - create samplesheet; assume '--name-mode noLaneSplit' by default
# location of production demultiplexing for deployment
FASTQDIR:=
# sequencing run name for deployment
RUNID:=
# location of production deployment for analysis
# NYUMC Big Purple directory location
PRODDIR:=/gpfs/data/molecpathlab/production/NGS580
deploy:
	@$(MAKE) check-runid
	@$(MAKE) check-fastqdir
	@repo_dir="$${PWD}" && \
	output_dir="$(PRODDIR)/$(RUNID)" && \
	echo ">>> Setting up new repo in location: $${output_dir}" && \
	git clone --recursive "$${repo_dir}" "$${output_dir}" && \
	cd "$${output_dir}" && \
	echo ">>> Linking input directory: $(FASTQDIR)" && \
	ln -s "$(FASTQDIR)" input && \
	echo ">>> Creating input fastq sheet" && \
	python generate-samplesheets.py --name-mode noLaneSplit '$(FASTQDIR)' && \
	echo ">>> Creating config file..." && \
	$(MAKE) config CONFIG_OUTPUT="$${output_dir}/config.json" && \
	printf ">>> NGS580 analysis directory prepared: $${output_dir}\n"

CONFIG_INPUT:=.config.json
CONFIG_OUTPUT:=config.json
$(CONFIG_OUTPUT):
	@echo ">>> Creating $(CONFIG_OUTPUT)"
	@cp "$(CONFIG_INPUT)" "$(CONFIG_OUTPUT)"

SAMPLESHEET:=
config: $(CONFIG_OUTPUT)
	@[ -n "$(RUNID)" ] && echo ">>> Updating runID config" && python config.py --update "$(CONFIG_OUTPUT)" --runID "$(RUNID)" || :
	@[ -n "$(SAMPLESHEET)" ] && echo ">>> Updating samplesheet config" && python config.py --update "$(CONFIG_OUTPUT)" --samplesheet "$(SAMPLESHEET)" || :
	@[ -n "$(FASTQDIR)" ] && echo ">>> Updating fastqDirs config" && python config.py --update "$(CONFIG_OUTPUT)" --fastqDirs "$(FASTQDIR)" || :

# generate a samplesheet for the analysis
samplesheet: check-fastqdir samples.analysis.tsv
	echo ">>> Generating samplesheet for fastq directory $(FASTQDIR)" && \
	python generate-samplesheets.py '$(FASTQDIR)' && \
	$(MAKE) config SAMPLESHEET=samples.analysis.tsv

# add an extra fastq directory to the samplesheet
# update the analysis config for the new directory
samplesheet-add-fastqdir: check-fastqdir
	@echo ">>> Adding fastq directory $(FASTQDIR) to the samplesheet" && \
	python generate-samplesheets.py --append '$(FASTQDIR)' && \
	old_fastqdirs="$$(python -c 'import json; fastq_dirs = json.load(open("config.json")).get("fastqDirs", None); print(" " .join(fastq_dirs) if fastq_dirs else "" )')" && \
	echo ">>> Old fastq directories were: $${old_fastqdirs}" && \
	new_fastqdirs="$(FASTQDIR) $${old_fastqdirs}" && \
	echo ">>> New fastq directories will be: $${new_fastqdirs}" && \
	python config.py --update "$(CONFIG_OUTPUT)" --fastqDirs $${new_fastqdirs}

# ~~~~~ UPDATE THIS REPO ~~~~~ #
# commands for bringing this directory's pipeline up to date
update: pull update-submodules update-nextflow

# pull latest version of repo
pull: remote
	@echo ">>> Updating repo"
	@git pull

# remove local Nextflow install and re-build
update-nextflow:
	@if [ -f nextflow ]; then \
	echo ">>> Removing old Nextflow" && \
	rm -f nextflow && \
	echo ">>> Reinstalling Nextflow" && \
	$(MAKE) install ; \
	else $(MAKE) install ; fi

# pull the latest version of all submodules
update-submodules: remote
	@echo ">>> Updating git submodules"
	@git submodule update --recursive --remote --init

# update the repo remote for ssh
remote-ssh:
	@git remote set-url origin "$(REMOTE_ssh)"
remote:
	@echo ">>> Setting git remote origin to $(REMOTE_http)"
	@git remote set-url origin "$(REMOTE_http)"


# ~~~~~ RUN PIPELINE ~~~~~ #
RESUME:=-resume

# try to automatically determine which 'run' recipe to use based on hostname
run:
	@if grep -q 'phoenix' <<<'$(HOSTNAME)'; then \
	$(MAKE) run-phoenix; \
	elif grep -q 'bigpurple' <<<'$(HOSTNAME)'; then \
	$(MAKE) run-bigpurple ; \
	else echo ">>> ERROR: could not automatically determine 'run' recipe to use, please consult the Makefile"; exit 1 ; \
	fi ; \

# run on phoenix in current session
run-phoenix: install
	@module unload java && module load java/1.8 && \
	log_file="logs/nextflow.$(TIMESTAMP).stdout.log" ; \
	echo ">>> Running Nextflow with stdout log file: $${log_file}" ; \
	./nextflow run main.nf -profile phoenix $(RESUME) -with-dag flowchart.dot $(EP) | tee -a "$${log_file}" ; \
	echo ">>> Nextflow completed, stdout log file: $${log_file}"

# run on NYU Big Purple HPC in current session
Q:=cpu_medium
run-bigpurple: install
	@log_file="logs/nextflow.$(TIMESTAMP).stdout.log" ; \
	echo ">>> Running Nextflow with stdout log file: $${log_file}" ; \
	./nextflow run main.nf -profile bigPurple $(RESUME) -with-dag flowchart.dot --queue $(Q) $(EP) | tee -a "$${log_file}" ; \
	echo ">>> Nextflow completed, stdout log file: $${log_file}"

# run locally default settings
run-local: install
	./nextflow run main.nf -profile local $(RESUME) -with-dag flowchart.dot $(EP)

# run on Power server
run-power: install
	source /shared/miniconda2/bin/activate /shared/biobuilds-2017.11 && \
	./nextflow run main.nf -profile power $(RESUME) -with-dag flowchart.dot $(EP)

# submit the parent Nextflow process to phoenix HPC as a qsub job
# `make submit-phoenix EP='--runID 180316_NB501073_0036_AH3VFKBGX5'`
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
	./nextflow run main.nf -profile phoenix $(RESUME) -with-dag flowchart.dot $(EP) & \
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
	@[ -z "$(REMOTE)" ] && printf "invalid REMOTE server specified: $(REMOTE)\n" && exit 1 || :
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

# compile flow chart
flowchart:
	[ -f flowchart.dot ] && dot flowchart.dot -Tpng -o flowchart.png || echo "file flowchart.dot not present"

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
#  run with `make finalize -j8`

# Nextflow "publishDir" directory of files to keep
publishDir:=output
# Nextflow "work" directory
workDir:=work
TRACEFILE:=trace.txt

# remove extraneous work dirs
# resolve publishDir output symlinks
# write work 'ls' files
# create work dir file stubs
finalize:
	$(MAKE) finalize-work-rm
	$(MAKE) finalize-output
	$(MAKE) finalize-work-ls
	$(MAKE) finalize-work-stubs

## ~~~ convert all symlinks to their linked items ~~~ ##
# symlinks in the publishDir to convert to files
publishDirLinks:=
FIND_publishDirLinks:=
ifneq ($(FIND_publishDirLinks),)
publishDirLinks:=$(shell find $(publishDir)/ -type l)
endif
finalize-output:
	@echo ">>> Converting symlinks in output dir '$(publishDir)' to their targets..."
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


## ~~~ write list of files in each subdir to file '.ls.txt' ~~~ ##
# subdirs in the 'work' dir
NXFWORKSUBDIRS:=
FIND_NXFWORKSUBDIRS:=
ifneq ($(FIND_NXFWORKSUBDIRS),)
NXFWORKSUBDIRS:=$(shell find "$(workDir)/" -maxdepth 2 -mindepth 2)
endif
# file to write 'ls' contents of 'work' subdirs to
LSFILE:=.ls.txt
finalize-work-ls:
	@echo ">>> Writing list of directory contents for each subdir in Nextflow work directory '$(workDir)'..."
	$(MAKE) finalize-work-ls-recurse FIND_NXFWORKSUBDIRS=1
finalize-work-ls-recurse: $(NXFWORKSUBDIRS)
# print the 'ls' contents of each subdir to a file, or delete the subdir
$(NXFWORKSUBDIRS):
	@ls_file="$@/$(LSFILE)" ; \
	echo ">>> Writing file list: $${ls_file}" ; \
	ls -1 "$@" > "$${ls_file}"
.PHONY: $(NXFWORKSUBDIRS)



## ~~~ replace all files in 'work' dirs with empty file stubs ~~~ ##
NXFWORKFILES:=
FIND_NXFWORKFILES:=
# files in work subdirs to keep
LSFILEREGEX:=\.ls\.txt
NXFWORKFILES:='.command.begin|.command.err|.command.log|.command.out|.command.run|.command.sh|.command.stub|.command.trace|.exitcode|$(LSFILE)'
NXFWORKFILESREGEX:='.*\.command\.begin\|.*\.command\.err\|.*\.command\.log\|.*\.command\.out\|.*\.command\.run\|.*\.command\.sh\|.*\.command\.stub\|.*\.command\.trace\|.*\.exitcode\|.*$(LSFILEREGEX)'
ifneq ($(FIND_NXFWORKFILES),)
NXFWORKFILES:=$(shell find -P "$(workDir)/" -type f ! -regex $(NXFWORKFILESREGEX))
endif
finalize-work-stubs:
	$(MAKE) finalize-work-stubs-recurse FIND_NXFWORKFILES=1
finalize-work-stubs-recurse: $(NXFWORKFILES)
$(NXFWORKFILES):
	@printf '>>> Creating file stub: $@\n' && rm -f "$@" && touch "$@"
.PHONY: $(NXFWORKFILES)


## ~~~ remove 'work' subdirs that are not in the latest trace file (e.g. most previous run) ~~~ ##
# subdirs in the 'work' dir
NXFWORKSUBDIRSRM:=
FIND_NXFWORKSUBDIRSRM:=
# regex from the hashes of tasks in the tracefile to match against work subdirs
HASHPATTERN:=
ifneq ($(FIND_NXFWORKSUBDIRSRM),)
NXFWORKSUBDIRSRM:=$(shell find "$(workDir)/" -maxdepth 2 -mindepth 2)
HASHPATTERN:=$(shell python -c 'import csv; reader = csv.DictReader(open("$(TRACEFILE)"), delimiter = "\t"); print("|".join([row["hash"] for row in reader]))')
endif
finalize-work-rm:
	@echo ">>> Removing subdirs in Nextflow work directory '$(workDir)' which are not included in Nextflow trace file '$(TRACEFILE)'..."
	$(MAKE) finalize-work-rm-recurse FIND_NXFWORKSUBDIRSRM=1
finalize-work-rm-recurse: $(NXFWORKSUBDIRSRM)
# remove the subdir if its not listed in the trace hashes
$(NXFWORKSUBDIRSRM):
	@if [ ! "$$(echo '$@' | grep -q -E "$(HASHPATTERN)"; echo $$? )" -eq 0 ]; then \
	echo ">>> Removing subdir: $@" ; \
	rm -rf "$@" ; \
	fi
.PHONY: $(NXFWORKSUBDIRSRM)
