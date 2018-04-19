SHELL:=/bin/bash
VAR=
none:

.PHONY: $(VAR)

check-Docker:
	docker --version > /dev/null 2>&1 || { echo "ERROR: 'docker' not found" && exit 1 ; }

build-all-Docker: check-Docker
	for item in $$(find . -maxdepth 2 -mindepth 2 -type f -name "Dockerfile" -exec dirname {} \; | cut -d '/' -f2- | xargs); \
	do \
	make -f Docker.makefile build VAR=$$item ; \
	done

pull-all-Docker: check-Docker
	for item in $$(find . -maxdepth 2 -mindepth 2 -type f -name "Dockerfile" -exec dirname {} \; | cut -d '/' -f2- | xargs); \
	do \
	make -f Docker.makefile pull VAR=$$item ; \
	done

convert-all-Docker: check-Docker
	for item in $$(find . -maxdepth 2 -mindepth 2 -type f -name "Dockerfile" -exec dirname {} \; | cut -d '/' -f2- | xargs); \
	do \
	make -f Singularity.makefile convert-Docker VAR=$$item ; \
	done

build-all-Singularity:
	for item in $$(find . -maxdepth 2 -mindepth 2 -type f -name "Singularity.*" -exec dirname {} \; | cut -d '/' -f2- | xargs); \
	do \
	make -f Singularity.makefile build VAR=$$item ; \
	done
