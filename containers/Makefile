# building Docker cotainers
none:

.PHONY: base fastqc trimmomatic sambamba bwa bedtools

base:
	cd base && docker build -t stevekm/ngs580-nf:base .

gatk:
	docker pull broadinstitute/gatk:4.0.1.2

bwa: base
	cd bwa && docker build -t stevekm/ngs580-nf:bwa .

fastqc: base
	cd fastqc && docker build -t stevekm/ngs580-nf:fastqc .

trimmomatic: base
	cd trimmomatic && docker build -t stevekm/ngs580-nf:trimmomatic .

sambamba: base
	cd sambamba && docker build -t stevekm/ngs580-nf:sambamba .

bedtools: base
	cd bedtools && docker build -t stevekm/ngs580-nf:bedtools .

variant-calling: gatk
	/bin/cp ../bin/annotate_vcf.sh variant-calling/
	cd variant-calling && docker build -t stevekm/ngs580-nf:variant-calling .

multiqc: base
	/bin/cp ../bin/multiqc.requirements.txt multiqc
	cd multiqc && docker build -t stevekm/ngs580-nf:multiqc .



build: base gatk bwa fastqc trimmomatic sambamba bedtools variant-calling multiqc

pull:
	docker pull broadinstitute/gatk:4.0.1.2
	docker pull stevekm/ngs580-nf:base
	docker pull stevekm/ngs580-nf:bwa
	docker pull stevekm/ngs580-nf:fastqc
	docker pull stevekm/ngs580-nf:trimmomatic
	docker pull stevekm/ngs580-nf:sambamba
	docker pull stevekm/ngs580-nf:bedtools
	docker pull stevekm/ngs580-nf:variant-calling
	docker pull stevekm/ngs580-nf:multiqc



test-base: base
	docker run --rm -ti stevekm/ngs580-base bash

test-bwa: bwa
	docker run --rm -ti stevekm/ngs580-nf:bwa bash

test-fastqc: fastqc
	docker run --rm -ti stevekm/ngs580-nf:fastqc bash

test-gatk: gatk
	docker run --rm -ti broadinstitute/gatk bash

test-trimmomatic: trimmomatic
	docker run --rm -ti stevekm/ngs580-nf:trimmomatic bash

test-sambamba: sambamba
	docker run --rm -ti stevekm/ngs580-nf:sambamba bash

test-bedtools: bedtools
	docker run --rm -ti stevekm/ngs580-nf:bedtools bash

test-variant-calling: variant-calling
	docker run --rm -ti stevekm/ngs580-nf:variant-calling bash

test-multiqc: multiqc
	docker run --rm -ti stevekm/ngs580-multiqc bash
