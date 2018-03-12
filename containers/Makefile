none:

base:
	cd ngs580-base && docker build -t stevekm/ngs580-base .

main: base
	cd ngs580-main && docker build -t stevekm/ngs580-main .

gatk:
	docker pull broadinstitute/gatk

bwa:
	cd ngs580-bwa && docker build -t stevekm/ngs580-bwa .

fastqc:
	cd ngs580-FastQC && docker build -t stevekm/ngs580-fastqc .

trimmomatic:
	cd ngs580-Trimmomatic && docker build -t stevekm/ngs580-trimmomatic .



build: main gatk bwa

test: base
	docker run --rm -ti stevekm/ngs580-base bash

test-main: main
	docker run --rm -ti stevekm/ngs580-main bash

test-bwa: bwa
	docker run --rm -ti stevekm/ngs580-bwa bash

test-fastqc: fastqc
	docker run --rm -ti stevekm/ngs580-fastqc bash

test-gatk: gatk
	docker run --rm -ti broadinstitute/gatk bash

test-trimmomatic: trimmomatic
	docker run --rm -ti stevekm/ngs580-trimmomatic bash
