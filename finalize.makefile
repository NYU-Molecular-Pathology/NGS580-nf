# run with 'make -f finalize.makefile -j4'
publishDir:=output
publishDirLinks:=$(shell find $(publishDir)/ -type l)
.PHONY: output $(publishDirLinks)

none:

output: $(publishDirLinks)

# replace symlinks in publishDir with the original file or dir
$(publishDirLinks):
	@ { \
	destination="$@"; \
	sourcepath="$$(python -c "import os; print(os.path.realpath('$@'))")" ; \
	if [ ! -e "$${sourcepath}" ]; then echo "ERROR: Source does not exist: $${sourcepath}"; \
	elif [ -f "$${sourcepath}" ]; then rsync -va --remove-source-files "$$sourcepath" "$$destination" ; \
	elif [ -d "$${sourcepath}" ]; then { \
	timestamp="$$(date +%s)" ; \
	tmpdir="$${destination}.$${timestamp}" ; \
	rsync -va --remove-source-files "$${sourcepath}/" "$${tmpdir}" && \
	rm -f "$${destination}" && \
	mv "$${tmpdir}" "$${destination}" ; } ; \
	fi ; }
