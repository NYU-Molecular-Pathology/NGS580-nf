# Makefile to run the pipeline
SHELL:=/bin/bash
export NXF_VER:=19.07.0
export NXF_ANSI_LOG=false
# extra params to pass for Nextflow in some recipes
EP:=
TIMESTAMP:=$(shell date +%s)
TIMESTAMP_str:=$(shell date +"%Y-%m-%d-%H-%M-%S")
ABSDIR:=$(shell python -c 'import os; print(os.path.realpath("."))')
DIRNAME:=$(shell python -c 'import os; print(os.path.basename(os.path.realpath(".")))')
REMOTE_ssh:=git@github.com:NYU-Molecular-Pathology/NGS580-nf.git
REMOTE_http:=https://github.com/NYU-Molecular-Pathology/NGS580-nf.git
.PHONY: annovar_db ref

# no default action
none:

# ~~~~~ SETUP PIPELINE ~~~~~ #
HOSTNAME:=$(shell echo $$HOSTNAME)
USER_HOME=$(shell echo "$$HOME")
USER_DATE:=$(shell date +%s)

# relative locations
# Nextflow "publishDir" directory of files to keep
publishDir:=output
# Nextflow "work" directory
workDir:=work
TRACEFILE:=trace.txt

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
	$(MAKE) config RUNID=demo FASTQDIR=NGS580-demo-data/tiny/fastq/ && \
	$(MAKE) samplesheet && \
	mv targets.bed "targets.bed.$(TIMESTAMP)" && \
	/bin/cp NGS580-demo-data/tiny/targets.bed . && \
	./update-samplesheets.py --tumor-normal-sheet NGS580-demo-data/tiny/samples.pairs.csv \
	--pairs-tumor-colname '#SAMPLE-T' \
	--pairs-normal-colname '#SAMPLE-N'
	# ./generate-samplesheets.py NGS580-demo-data/tiny/fastq/ && \
	# mv targets.bed "targets.bed.$(TIMESTAMP)" && \
	# /bin/cp NGS580-demo-data/tiny/targets.bed .
	# ./update-samplesheets.py --tumor-normal-sheet NGS580-demo-data/tiny/samples.pairs.csv \
	# --pairs-tumor-colname '#SAMPLE-T' \
	# --pairs-normal-colname '#SAMPLE-N'
	#
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
FASTQDIRS:=
# samplesheet used for sequencing run demultiplexing
DEMUX_SAMPLESHEET:=
DEMUX_SAMPLESHEET_output:=demux-samplesheet.csv
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
	[ -e "$(DEMUX_SAMPLESHEET)" ] && /bin/cp "$(DEMUX_SAMPLESHEET)" "$(DEMUX_SAMPLESHEET_output)" || : \
	echo ">>> Creating input fastq sheet" && \
	python generate-samplesheets.py --name-mode "$(NAMEMODE)" '$(FASTQDIR)' && \
	echo ">>> Creating config file..." && \
	$(MAKE) config CONFIG_OUTPUT="$${output_dir}/config.json" && \
	[ -e "$(DEMUX_SAMPLESHEET_output)" ] && echo ">>> Adding demux samplesheet to config" && $(MAKE) config DEMUX_SAMPLESHEET="$(DEMUX_SAMPLESHEET_output)" CONFIG_OUTPUT="$${output_dir}/config.json" || : \
	[ -e "$(DEMUX_SAMPLESHEET_output)" ] && $(MAKE) pairs PAIRS_SHEET="$(DEMUX_SAMPLESHEET_output)" PAIRS_MODE=demux
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
	@[ -n "$(DEMUX_SAMPLESHEET)" ] && echo ">>> Updating demultiplexing samplesheet config" && python config.py --update "$(CONFIG_OUTPUT)" --demux-samplesheet "$(DEMUX_SAMPLESHEET)" || :
	@[ -n "$(FASTQDIR)" ] && echo ">>> Updating fastqDirs config" && python config.py --update "$(CONFIG_OUTPUT)" --fastqDirs "$(FASTQDIR)" || :
	@[ -n "$(FASTQDIRS)" ] && echo ">>> Adding fastq dirs to config" && python config.py --update "$(CONFIG_OUTPUT)" --fastqDirs $(FASTQDIRS) || :

config-add-fastqdirs:
	@if [ ! -z "$(FASTQDIRS)" ]; then \
	for fastqdir in $(FASTQDIRS); do \
	echo ">>> Adding $$fastqdir to $(CONFIG_OUTPUT)" ; \
	$(MAKE) config FASTQDIR="$$fastqdir" ; \
	done; else \
	echo ">>> ERROR: no FASTQDIRS passed" ; \
	fi

# generate a samplesheet for the analysis based on the configs
NAMEMODE:=noLaneSplit
SAMPLESHEET_OUTPUT:=samples.analysis.tsv
samplesheet:
	@echo ">>> Getting fastqdirs from config file: $(CONFIG_OUTPUT)" && \
	fastqdirs="$$(python -c 'import json; fastq_dirs = json.load(open("$(CONFIG_OUTPUT)")).get("fastqDirs", None); print(" " .join(fastq_dirs) if fastq_dirs else "" )')" && \
	echo ">>> loaded fastqdirs: $${fastqdirs}" && \
	echo ">>> Generating samplesheet '$(SAMPLESHEET_OUTPUT)' for fastqdirs" && \
	python generate-samplesheets.py $(EP) --samples-analysis-tsv "$(SAMPLESHEET_OUTPUT)" --name-mode "$(NAMEMODE)" $${fastqdirs}

# update the samplesheet with the sample pairs information
# requires samplesheet to exist
# default to 'sns' style sample pairs .csv file
PAIRS_SHEET:=samples.pairs.csv
# type of sheet to parse pairs from;
# 'sns': Igor sns pipeline format
# 'demux': Illumina demultiplexing samplesheet with extra 'Paired_Normal' column added to [Data] section
PAIRS_MODE:=sns
pairs:
	@if [ ! -e "$(SAMPLESHEET_OUTPUT)" ]; then $(MAKE) samplesheet; fi && \
	if [ ! -e "$(PAIRS_SHEET)" ]; then echo ">>> ERROR: PAIRS_SHEET does not exist: $(PAIRS_SHEET)"; exit 1; fi && \
	if [ "$(PAIRS_MODE)" == "sns" ]; then \
	echo ">>> Updating samplesheet with sample pairs from sns style sheet" ; \
	python update-samplesheets.py --tumor-normal-sheet "$(PAIRS_SHEET)" --pairs-tumor-colname '#SAMPLE-T' --pairs-normal-colname '#SAMPLE-N' ; \
	elif [ "$(PAIRS_MODE)" == "demux" ]; then \
	echo ">>> Updating samplesheet with sample pairs from demultiplexing style sheet" ; \
	python bin/demux2tumor_normal_sheet.py "$(PAIRS_SHEET)" samples.tumor.normal.csv && \
	python update-samplesheets.py --tumor-normal-sheet samples.tumor.normal.csv ; \
	else echo ">>> ERROR: PAIRS_MODE not recognized: $(PAIRS_MODE)"; exit 1; fi

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
# enable 'resume' functionality; set this to '' to disable
RESUME:=-resume
LOGDIR:=logs
LOGDIRABS:=$(shell python -c 'import os; print(os.path.realpath("$(LOGDIR)"))')
LOGID:=$(TIMESTAMP)
LOGFILEBASE:=log.$(LOGID).out
LOGFILE:=$(LOGDIR)/$(LOGFILEBASE)
# start run with logging
run: backup
	@log_file="$(LOGDIR)/nextflow.$(LOGID).stdout.log" ; \
	echo ">>> Running with stdout log file: $(LOGFILE)" ; \
	$(MAKE) run-recurse 2>&1 | tee -a "$(LOGFILE)" ; \
	echo ">>> Run completed at $$(date +"%Y-%m-%d %H:%M:%S"), stdout log file: $(LOGFILE)"

# try to automatically determine which 'run' recipe to use based on hostname
run-recurse:
	@if grep -q 'phoenix' <<<'$(HOSTNAME)'; then echo  ">>> Running run-phoenix"; $(MAKE) run-phoenix ; \
	elif grep -q 'bigpurple' <<<'$(HOSTNAME)'; then echo ">>> Running run-bigpurple"; $(MAKE) run-bigpurple ; \
	else echo ">>> ERROR: could not automatically determine 'run' recipe to use, please consult the Makefile"; exit 1 ; fi ;

# run on phoenix in current session
run-phoenix: install
	@module unload java && module load java/1.8 && \
	./nextflow run main.nf -profile phoenix $(RESUME) -with-dag flowchart.dot $(EP)

# run on NYU Big Purple HPC
SLURM_recurse:=
# HPC queue/partition to submit to by default
Q:=cpu_short
# total number of CPUs on Intellispace queue
Intellispace_queue_max:=80
Intellispace_queue_name:=intellispace
ifneq ($(SLURM_recurse),)
# check the count of intellispace queue cpus used/requested
Intellispace_queue_cpus:=$(shell squeue -p intellispace -o "%C" | tail -n +2 | paste -sd+ | bc)
# try to automatically detect a SLURM queue with idle nodes to submit to
AutoQ_SLURM_idle:=$(shell sinfo -N -O nodelist,partition,statelong | grep 'idle' | grep -v 'data_mover' | grep -v 'dev' | grep -v 'fn_' | grep -v 'cpu_long' | tr -s '[:space:]' | cut -d ' ' -f2 | sort -u | head -1)
# detect which 'mixed' queue has the most open nodes
AutoQ_SLURM_mixed:=$(shell sinfo -N -O nodelist,partition,statelong | grep 'mixed' | grep -v 'data_mover' | grep -v 'dev' | grep -v 'fn_' | grep -v 'cpu_long' | tr -s '[:space:]' | cut -d ' ' -f2 | sort | uniq -c | sort -k 1nr | head -1 | tr -s '[:space:]' | cut -d ' ' -f3)
# first option: use queue with most idle nodes; empty if none exist
ifneq ($(AutoQ_SLURM_idle),)
Q:=$(AutoQ_SLURM_idle)
# second option: use queue with most mixed nodes; empty if none exist
else ifneq ($(AutoQ_SLURM_mixed),)
Q:=$(AutoQ_SLURM_mixed)
# third option: use intellispace queue if it has less than the max number of CPUs used/requested already
else ifneq ($(shell test $(Intellispace_queue_cpus) -lt $(Intellispace_queue_max); echo $$?),0)
Q:=$(Intellispace_queue_name)
endif
endif
test-q:
	@$(MAKE) test-q-recurse SLURM_recurse=1
test-q-recurse:
	@echo "Q: $(Q), AutoQ_SLURM_idle: $(AutoQ_SLURM_idle), AutoQ_SLURM_mixed: $(AutoQ_SLURM_mixed), Intellispace_queue_cpus: $(Intellispace_queue_cpus)"
run-bigpurple:
	$(MAKE) run-bigpurple-recurse SLURM_recurse=1
run-bigpurple-recurse: Q_JSON:=/gpfs/home/kellys04/molecpathlab/pipelines/queue-stats/slurm.json
run-bigpurple-recurse: export NXF_DEBUG=3
run-bigpurple-recurse: install
	./nextflow -trace nextflow.executor run main.nf -profile bigPurple $(RESUME) -with-dag flowchart.dot --queue_json "$(Q_JSON)" $(EP)
	$(MAKE) fix-permissions fix-group
# --queue "$(Q)" # try using the queue JSON instead

# run locally default settings
run-local: install
	./nextflow run main.nf -profile local $(RESUME) -with-dag flowchart.dot $(EP)

# run on Power server
run-power: install
	source /shared/miniconda2/bin/activate /shared/biobuilds-2017.11 && \
	./nextflow run main.nf -profile power $(RESUME) -with-dag flowchart.dot $(EP)


# submit the parent Nextflow process to phoenix HPC as a cluster job
SUBJOBNAME:=NGS580-$(DIRNAME)
SUBLOG:=$(LOGDIRABS)/slurm-%j.$(LOGFILEBASE)
SUBQ:=intellispace
SUBTIME:=--time=5-00:00:00
SUBTHREADS:=8
SUBMEM:=48G
SUBEP:=
NXF_NODEFILE:=.nextflow.node
NXF_JOBFILE:=.nextflow.jobid
NXF_PIDFILE:=.nextflow.pid
NXF_SUBMIT:=.nextflow.submitted
NXF_SUBMITLOG:=.nextflow.submitted.log
REMOTE:=
PID:=
SUBSCRIPT:=submit.bigpurple.sh
# check for an HPC submission lock file, then try to determine the submission recipe to use
submit:
	@if [ -e "$(NXF_SUBMIT)" ]; then echo ">>> ERROR: Workflow locked by $(NXF_SUBMIT); has an instance of the pipeline has already been submitted?"; exit 1 ; \
	else \
	if grep -q 'phoenix' <<<'$(HOSTNAME)'; then echo  ">>> Submission for phoenix not yet configured";  \
	elif grep -q 'bigpurple' <<<'$(HOSTNAME)'; then echo ">>> Running submit-bigpurple"; $(MAKE) submit-bigpurple ; \
	else echo ">>> ERROR: could not automatically determine 'submit' recipe to use, please consult the Makefile"; exit 1 ; fi ; \
	fi

# generate a SLURM sbatch script to start the pipeline in a SLURM job
# sets a submission lock file
# NOTE: Nextflow locks itself from concurrent instances but need to lock against multiple 'make submit'
# catches 'scancel' commands and propagates them up to Nextflow for clean shutdown
$(SUBSCRIPT):
	@printf '' > $(SUBSCRIPT)
	@chmod +x $(SUBSCRIPT)
	@echo '#!/bin/bash -x' >> $(SUBSCRIPT)
	@echo '#SBATCH -D $(ABSDIR)' >> $(SUBSCRIPT)
	@echo '#SBATCH -o $(SUBLOG)' >> $(SUBSCRIPT)
	@echo '#SBATCH -J $(SUBJOBNAME)' >> $(SUBSCRIPT)
	@echo '#SBATCH -p $(SUBQ)' >> $(SUBSCRIPT)
	@echo '#SBATCH $(SUBTIME)' >> $(SUBSCRIPT)
	@echo '#SBATCH --ntasks-per-node=1' >> $(SUBSCRIPT)
	@echo '#SBATCH -c $(SUBTHREADS)' >> $(SUBSCRIPT)
	@echo '#SBATCH --mem $(SUBMEM)' >> $(SUBSCRIPT)
	@echo '#SBATCH --export=HOSTNAME' >> $(SUBSCRIPT)
	@echo 'touch $(NXF_SUBMIT)' >> $(SUBSCRIPT)
	@echo 'get_pid(){ head -1 $(NXF_PIDFILE); }' >> $(SUBSCRIPT)
	@echo 'rm_submit(){ echo ">>> trap: rm_submit" ; [ -e $(NXF_SUBMIT) ] && rm -f $(NXF_SUBMIT) || : ; }' >> $(SUBSCRIPT)
	@echo 'wait_pid(){ local pid=$$1 ; while kill -0 $$pid; do echo waiting for process $$pid to end ; sleep 1 ; done ; }' >> $(SUBSCRIPT)
	@echo 'nxf_kill(){ rm_submit ; echo ">>> trap: nxf_kill" && pid=$$(get_pid) && kill $$pid && wait_pid $$pid ; }' >> $(SUBSCRIPT)
	@echo 'trap nxf_kill HUP' >> $(SUBSCRIPT)
	@echo 'trap nxf_kill INT' >> $(SUBSCRIPT)
	@echo 'trap nxf_kill EXIT' >> $(SUBSCRIPT)
	@echo 'make submit-bigpurple-run TIMESTAMP=$(TIMESTAMP) $(SUBEP)' >> $(SUBSCRIPT)
.PHONY: $(SUBSCRIPT)

# submit on Big Purple using SLURM
submit-bigpurple: $(SUBSCRIPT)
	@sbatch $(SUBSCRIPT) | tee >(sed 's|[^[:digit:]]*\([[:digit:]]*\).*|\1|' > '$(NXF_JOBFILE)')

# run inside a SLURM sbatch
# store old pid and node entries in a backup file in case things get messy
# need to manually set the HOSTNAME here because it changes inside SLURM job
# TODO: come up with a better method for this ^^
submit-bigpurple-run:
	if [ -e "$(NXF_NODEFILE)" -a -e "$(NXF_PIDFILE)" ]; then paste "$(NXF_NODEFILE)" "$(NXF_PIDFILE)" >> $(NXF_SUBMITLOG); fi ; \
	echo "$${SLURMD_NODENAME}" > "$(NXF_NODEFILE)" && \
	$(MAKE) run HOSTNAME="bigpurple" LOGID="$(TIMESTAMP)"

# issue an interupt signal to a process running on a remote server
# e.g. Nextflow running in a qsub job on a compute node
kill: PID=$(shell head -1 "$(NXF_PIDFILE)")
kill: REMOTE=$(shell head -1 "$(NXF_NODEFILE)")
kill: $(NXF_NODEFILE) $(NXF_PIDFILE)
	ssh "$(REMOTE)" 'kill $(PID)'

# location to files for making HapMap pool
HAPMAP_POOL_SHEET:=samples.hapmap.tsv
hapmap-pool: $(HAPMAP_POOL_SHEET)
	./nextflow run hapmap-pool.nf -profile hapmap_pool $(RESUME)

# location for files for making CNV Pool
CNV_POOL_SHEET:=samples.cnv.tsv
cnv-pool: $(CNV_POOL_SHEET)
	./nextflow run cnv-pool.nf -profile cnv_pool $(RESUME)

# save a record of the most recent Nextflow run completion
PRE:=
RECDIR:=recorded-runs/$(PRE)$(DIRNAME)_$(TIMESTAMP_str)
STDOUTLOGPATH:=
STDOUTLOG:=
ALL_LOGS:=
_RECORD:=
# task in the workdir e.g. '64/c111'
TASK:=
TASKDIR:=
ifneq ($(_RECORD),)
STDOUTLOGPATH:=$(shell ls -d -1t $(LOGDIR)/log.*.out | head -1 | python -c 'import sys, os; print(os.path.realpath(sys.stdin.readlines()[0].strip()))' )
STDOUTLOG:=$(shell basename "$(STDOUTLOGPATH)")
ALL_LOGS:=$(shell find "$(LOGDIR)" -type f -name '*$(STDOUTLOG)*')
ifneq ($(TASK),)
TASKDIR:=$(shell find "$(workDir)/" -mindepth 2 -maxdepth 2 -type d -path "*$(TASK)*" | head -1)
endif
endif
record:
	$(MAKE) record-recurse _RECORD=1
record-recurse:
	@echo ">>> Recording logs to $(RECDIR)" && \
	mkdir -p "$(RECDIR)" && \
	cp -a *.html trace.txt .nextflow.log main.nf nextflow.config config.json "$(RECDIR)/" ; \
	for item in $(ALL_LOGS); do cp -a "$${item}" "$(RECDIR)/"; done ; \
	if [ ! -z "$(TASKDIR)" -a -d "$(TASKDIR)" ]; then \
	echo ">>> Copying task dir: $(TASKDIR)" && \
	mkdir -p "$(RECDIR)/$(TASKDIR)" && \
	cp -a "$(TASKDIR)" "$(RECDIR)/$$(dirname $(TASKDIR))" ; \
	fi ; \
	echo ">>> Copied execution reports and logs to: $(RECDIR)"

# backup the analysis output
BACKUP_DIR:=$(publishDir)-backup/$(TIMESTAMP_str)
$(BACKUP_DIR):
	mkdir -p "$(BACKUP_DIR)"
backup:
	if [ -d "$(publishDir)" ] ; then \
	$(MAKE) $(BACKUP_DIR) BACKUP_DIR=$(BACKUP_DIR) ; \
	mv "$(publishDir)" "$(BACKUP_DIR)" ; fi


# fix group executable permissions
fix-perm:
	find . -type f -name "*.py" -o -name "*.R" ! -path "*/work/*" ! -path "*/output/*" -exec chmod -v g+x {} \;


# fix permissions on this directory
# make all executables group executable
# make all dirs full group accessible
# make all files group read/write
fix-permissions:
	@find . -type f -executable -exec chmod ug+X {} \;
	@find . -type d -exec chmod ug+rwxs {} \;
	@find . -type f -exec chmod ug+rw {} \;

USERGROUP:=molecpathlab
fix-group:
	@find . ! -group "$(USERGROUP)" -exec chgrp "$(USERGROUP)" {} \;





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
TRACE_PATTERN_FILE:=.trace.hash.txt
NXFWORKSUBDIRSRM:=
FIND_NXFWORKSUBDIRSRM:=
# regex from the hashes of tasks in the tracefile to match against work subdirs
HASHPATTERN:=
ifneq ($(FIND_NXFWORKSUBDIRSRM),)
NXFWORKSUBDIRSRM:=$(shell find "$(workDir)/" -maxdepth 2 -mindepth 2 )
# HASHPATTERN:=$(shell python -c 'import csv; reader = csv.DictReader(open("$(TRACEFILE)"), delimiter = "\t"); print("|".join([row["hash"] for row in reader]))')
endif
finalize-work-rm:
	@echo ">>> Removing subdirs in Nextflow work directory '$(workDir)' which are not included in Nextflow trace file '$(TRACEFILE)'..."
	$(MAKE) $(TRACE_PATTERN_FILE) finalize-work-rm-recurse FIND_NXFWORKSUBDIRSRM=1

# need to write out hashes to file because it gets too long to use as CLI arg
$(TRACE_PATTERN_FILE):
	@echo ">>> Making trace hash file: $(TRACE_PATTERN_FILE)"
	@python -c 'import csv; \
	reader = csv.DictReader(open("$(TRACEFILE)"), delimiter = "\t"); \
	fout = open("$(TRACE_PATTERN_FILE)", "w"); \
	fout.write("\n".join([row["hash"] for row in reader])) ; fout.close(); \
	'

finalize-work-rm-recurse: $(NXFWORKSUBDIRSRM)

$(NXFWORKSUBDIRSRM):
	@pattern="$$(echo $(@) | sed -e 's|$(workDir)/||g' | cut -c 1-9)" ; \
	if [ ! "$$( grep -q $${pattern} '$(TRACE_PATTERN_FILE)'; echo $$?)" -eq 0 ]; then \
	echo ">>> Removing subdir: $@" ; \
	rm -rf "$@" ; \
	fi

# remove the subdir if its not listed in the trace hashes
# $(NXFWORKSUBDIRSRM):
# 	@if [ ! "$$(echo '$@' | grep -q -E "$(HASHPATTERN)"; echo $$? )" -eq 0 ]; then \
# 	echo ">>> Removing subdir: $@" ; \
# 	rm -rf "$@" ; \
# 	fi
.PHONY: $(NXFWORKSUBDIRSRM) $(TRACE_PATTERN_FILE)
