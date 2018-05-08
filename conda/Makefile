SHELL:=/bin/bash
CONDA:=/shared/miniconda3

all: standard custom

standard:
	@for item in $$(find . -name env.yml -exec dirname {} \; | cut -d '/' -f2- | xargs); do \
	make -f conda.makefile create VAR="$$item" CONDA="$(CONDA)" ; \
	done

custom: 
	@for item in r-3.4.2 variant-calling-0.0.2; do \
	make -f conda.makefile create VAR="$$item" CONDA="$(CONDA)" ; \
	done

