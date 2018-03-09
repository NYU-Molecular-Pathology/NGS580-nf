none:

base:
	cd ngs580-base && docker build -t stevekm/ngs580-base .

main: base
	cd ngs580-main && docker build -t stevekm/ngs580-main .

gatk:
	docker pull broadinstitute/gatk

build: main gatk

test: base
	docker run --rm -ti stevekm/ngs580-base bash

test-main: main
	docker run --rm -ti stevekm/ngs580-main bash
