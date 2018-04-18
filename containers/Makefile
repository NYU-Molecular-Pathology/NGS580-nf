SHELL:=/bin/bash
VAR=
none:

.PHONY: $(VAR)

check-docker:
	docker --version > /dev/null 2>&1 || { echo "ERROR: 'docker' not found" && exit 1 ; }

build-all-Docker: check-docker
	for item in $$(find . -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | xargs); \
	do \
	make -f Docker.makefile build VAR=$$item ; \
	done

pull-all-Docker: check-docker
	for item in $$(find . -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | xargs); \
	do \
	make -f Docker.makefile pull VAR=$$item ; \
	done

convert-all-Docker: check-docker
	for item in $$(find . -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | xargs); \
	do \
	make -f Singularity.makefile convert-Docker VAR=$$item ; \
	done
