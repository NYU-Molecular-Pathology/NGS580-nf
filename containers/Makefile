none:

base:
	cd ngs580-base && docker build -t stevekm/ngs580-base .

main: base
	cd ngs580-main && docker build -t stevekm/ngs580-main .

gatk:
	docker pull broadinstitute/gatk

bwa: base
	cd ngs580-bwa && docker build -t stevekm/ngs580-bwa .

build: main gatk bwa

test: base
	docker run --rm -ti stevekm/ngs580-base bash

test-main: main
	docker run --rm -ti stevekm/ngs580-main bash

test-bwa: bwa
	docker run --rm -ti stevekm/ngs580-bwa bash

test-gatk: gatk
	docker run --rm -ti broadinstitute/gatk bash
