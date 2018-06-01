#  This Makfile contains commands to 'finalize' the Nextflow 'work' 
#  directory after the completion of a pipeline
#  Work dir contains subdirectories with intermediary files created when executing pipeline
#  along with command execution logs, etc
#  logs should be saved, along with a record of the contents of each subdir
#  files in the 'publishDir' are often symlinked to items in the 'work' dir;
#  these symlinks should be replaced with copies of the original file
#  all files in the 'work' dir should be removed, except for the log files
# 
#  Makefile configured for parallel processing of files
# 
#  run with `make -f finalize.makefile <recipe> -j8`
# 

# Nextflow "publishDir" directory of files to keep
publishDir:=output
publishDirLinks:=$(shell find $(publishDir)/ -type l)
# Nextflow "work" directory of items to be removed
workDir:=work
NXFWORKSUBDIRS:=$(shell find "$(workDir)/" -maxdepth 2 -mindepth 2) 
# file to write contents of 'work' subdirs to
LSFILE:=.ls.txt
LSFILEREGEX:=\.ls\.txt
# files in work subdirs to keep
NXFWORKFILES:='.command.begin|.command.err|.command.log|.command.out|.command.run|.command.sh|.command.stub|.command.trace|.exitcode|$(LSFILE)'
NXFWORKFILESREGEX:='.*\.command\.begin\|.*\.command\.err\|.*\.command\.log\|.*\.command\.out\|.*\.command\.run\|.*\.command\.sh\|.*\.command\.stub\|.*\.command\.trace\|.*\.exitcode\|.*$(LSFILEREGEX)'
# Nextflow 'trace' file with record of completed pipeline tasks
TRACEFILE:=trace-NGS580.txt
HASHPATTERN:=$(shell python -c 'import csv; reader = csv.DictReader(open("$(TRACEFILE)"), delimiter = "\t"); print("|".join([row["hash"] for row in reader]))')
.PHONY: output $(publishDirLinks) $(NXFWORKSUBDIRS)

none:

# replace symlinks in publishDir with the original item
output: $(publishDirLinks)
$(publishDirLinks):
	@ { \
	destination="$@"; \
	sourcepath="$$(python -c "import os; print(os.path.realpath('$@'))")" ; \
	if [ ! -e "$${sourcepath}" ]; then echo "ERROR: Source does not exist: $${sourcepath}"; \
	elif [ -f "$${sourcepath}" ]; then rsync -va "$$sourcepath" "$$destination" ; \
	elif [ -d "$${sourcepath}" ]; then { \
	timestamp="$$(date +%s)" ; \
	tmpdir="$${destination}.$${timestamp}" ; \
	rsync -va "$${sourcepath}/" "$${tmpdir}" && \
	rm -f "$${destination}" && \
	mv "$${tmpdir}" "$${destination}" ; } ; \
	fi ; }

# print the 'ls' contents of each subdir to a file, or delete the subdir 
ls: $(NXFWORKSUBDIRS)
$(NXFWORKSUBDIRS): 
	@if [ "$$(echo '$@' | grep -q -E "$(HASHPATTERN)"; echo $$? )" -eq 0 ]; then \
	ls_file="$@/$(LSFILE)" ; \
	ls -1 "$@" > "$${ls_file}" ; \
	else \
	echo "$@ is not in the trace file, deleting..." ; rm -rf "$@" ; \
	fi

# remove all files in the work dir TODO: figure out how to handle dirs
# test: NXFWORKOUTPUTREMOVE:=$(shell find "$(workDir)/" -mindepth 3 -maxdepth 3 -type f ! -regex $(NXFWORKFILESREGEX))
# # test: .PHONY: $(NXFWORKOUTPUTREMOVE)
# test: $(NXFWORKOUTPUTREMOVE)
# $(NXFWORKOUTPUTREMOVE):
# 	echo "do thing with $@ here"
# 	

# FIND_FILES :=
# FILES :=
# ifneq ($(FIND_FILES),)
# FILES := $(shell find "$(workDir)/" -mindepth 3 -maxdepth 3 -type f ! -regex $(NXFWORKFILESREGEX))
# endif
# FOO:=baz
# test:
# 	$(MAKE) test-recurse FIND_FILES=1 

# test-recurse: $(FILES)
# $(FILES):
# 	echo "do thing with file $@ here; $(FOO)"
# .PHONY: test test-recurse $(FILES)