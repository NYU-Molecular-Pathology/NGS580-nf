SHELL:=/bin/bash
VAR=
none:

.PHONY: $(VAR)

# ~~~~~ DOCKER ~~~~~ #
check-docker:
	docker --version > /dev/null 2>&1 || { echo "ERROR: 'docker' not found" && exit 1 ; }

build-all-docker: check-docker
	for item in $$(find . -maxdepth 1 -mindepth 1 -type d -exec basename {} \; | xargs); \
	do \
	echo "$$item"; \
	make -f Docker.makefile build VAR=$$item ; \
	done
