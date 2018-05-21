SHELL:=/bin/bash
CONDA:=/shared/miniconda3

all: standard custom

standard:
	@for item in $$(find . -name env.yml -exec dirname {} \; | cut -d '/' -f2- | xargs); do \
	if [ ! "$$( source "$(CONDA)/bin/activate" "$${item}" >/dev/null 2>&1 ; echo $$? )" -eq 0 ] ; then \
	make -f conda.makefile create VAR="$$item" CONDA="$(CONDA)" ; \
	else printf ">>> Skipping env that already exists: $${item}\n" ; \
	fi ; \
	done

custom: 
	@for item in r-3.4.2 variant-calling-0.0.2 annovar-150617; do \
	if [ ! "$$( source "$(CONDA)/bin/activate" "$${item}" >/dev/null 2>&1 ; echo $$? )" -eq 0 ] ; then \
	make -f conda.makefile $$item CONDA="$(CONDA)" ; \
	else printf ">>> Skipping env that already exists: $${item}\n" ; \
	fi ; \
	done

# run with 'make remove -j4'
ENVNAMES:=$(shell find . -name env.yml -exec dirname {} \; | cut -d '/' -f2- | xargs)
.PHONY: $(ENVNAMES) r-3.4.2 variant-calling-0.0.2 annovar-150617
remove: $(ENVNAMES) r-3.4.2 variant-calling-0.0.2 annovar-150617
$(ENVNAMES) r-3.4.2 variant-calling-0.0.2 annovar-150617:
	env="$@" ; \
	if [ "$$( source "$(CONDA)/bin/activate" "$${env}" >/dev/null 2>&1 ; echo $$? )" -eq 0 ] ; then \
	make -f conda.makefile remove-env VAR="$$env" CONDA="$(CONDA)" ; \
	else printf ">>> Skipping env that doesnt exist: $${env}\n" ; \
	fi ; \
	if [ -d "$(CONDA)/envs/$${env}" ]; then rm -rf "$(CONDA)/envs/$${env}"; fi
