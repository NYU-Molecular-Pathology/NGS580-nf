none:

.PHONY: base fastqc trimmomatic

base:
	cd base && docker build -t stevekm/ngs580-nf:base .

# deprecated ...
main: base
	cd ngs580-main && docker build -t stevekm/ngs580-main .

gatk:
	docker pull broadinstitute/gatk

bwa:
	cd ngs580-bwa && docker build -t stevekm/ngs580-bwa .

fastqc: base
	cd fastqc && docker build -t stevekm/ngs580-nf:fastqc .

trimmomatic: base
	cd trimmomatic && docker build -t stevekm/ngs580-nf:trimmomatic  .

sambamba: base
	cd ngs580-sambamba && docker build -t stevekm/ngs580-sambamba .

bedtools: base
	cd ngs580-bedtools && docker build -t stevekm/ngs580-bedtools .

variant-calling:
	/bin/cp ../bin/annotate_vcf.sh ngs580-variant-calling/
	cd ngs580-variant-calling && docker build -t stevekm/ngs580-variant-calling .

multiqc: base
	/bin/cp ../bin/multiqc.requirements.txt ngs580-multiqc
	cd ngs580-multiqc && docker build -t stevekm/ngs580-multiqc .



all: base gatk bwa fastqc trimmomatic sambamba bedtools variant-calling multiqc




test-base: base
	docker run --rm -ti stevekm/ngs580-base bash

test-main: main
	docker run --rm -ti stevekm/ngs580-main bash

test-bwa: bwa
	docker run --rm -ti stevekm/ngs580-bwa bash

test-fastqc: fastqc
	docker run --rm -ti stevekm/ngs580-nf:fastqc bash

test-gatk: gatk
	docker run --rm -ti broadinstitute/gatk bash

test-trimmomatic: trimmomatic
	docker run --rm -ti stevekm/ngs580-nf:trimmomatic bash

test-sambamba: sambamba
	docker run --rm -ti stevekm/ngs580-sambamba bash

test-bedtools: bedtools
	docker run --rm -ti stevekm/ngs580-bedtools bash

test-variant-calling: variant-calling
	docker run --rm -ti stevekm/ngs580-variant-calling bash

test-multiqc: multiqc
	docker run --rm -ti stevekm/ngs580-multiqc bash
