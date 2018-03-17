none:

.PHONY: base fastqc trimmomatic sambamba bwa bedtools msisensor R-3.3.0 delly2-0.7.7


# ~~~~~ BUILD DOCKER CONTAINERS ~~~~~ #
base:
	cd base && docker build -t stevekm/ngs580-nf:base .

gatk:
	docker pull broadinstitute/gatk:4.0.1.2

bwa: base
	cd bwa-0.7.17 && docker build -t stevekm/ngs580-nf:bwa-0.7.17 .

fastqc: base
	cd fastqc-0.11.7 && docker build -t stevekm/ngs580-nf:fastqc-0.11.7 .

trimmomatic: base
	cd trimmomatic-0.36 && docker build -t stevekm/ngs580-nf:trimmomatic-0.36 .

sambamba: base
	cd sambamba-0.6.6 && docker build -t stevekm/ngs580-nf:sambamba-0.6.6 .

bedtools: base
	cd bedtools-2.26.0 && docker build -t stevekm/ngs580-nf:bedtools-2.26.0 .

variant-calling1: gatk
	/bin/cp ../bin/annotate_vcf.sh variant-calling-0.0.1/
	cd variant-calling-0.0.1 && docker build -t stevekm/ngs580-nf:variant-calling-0.0.1 .

variant-calling2:
	/bin/cp ../bin/annotate_vcf.sh variant-calling-0.0.2/
	cd variant-calling-0.0.2 && docker build -t stevekm/ngs580-nf:variant-calling-0.0.2 .

variant-calling: variant-calling2

multiqc: base
	/bin/cp ../bin/multiqc.requirements.txt multiqc-1.4/
	cd multiqc-1.4 && docker build -t stevekm/ngs580-nf:multiqc-1.4 .

msisensor: base
	cd msisensor-0.2 && docker build -t stevekm/ngs580-nf:msisensor-0.2 .

R-3.2.3: base
	/bin/cp ../bin/install.R R-3.2.3/
	cd R-3.2.3 && docker build -t stevekm/ngs580-nf:R-3.2.3 .

R-3.3.2: base
	/bin/cp ../bin/install.R R-3.3.2/
	cd R-3.3.2 && docker build -t stevekm/ngs580-nf:R-3.3.2 .

R: R-3.3.2

delly2-0.7.7: base
	cd delly2-0.7.7 && docker build -t stevekm/ngs580-nf:delly2-0.7.7 .

delly2: delly2-0.7.7

annovar-150617: base
	/bin/cp ../bin/annotate_vcf.sh annovar-150617/
	cd annovar-150617 && docker build -t stevekm/ngs580-nf:annovar-150617 .

annovar: annovar-150617

# ~~~~~~~ SETUP DOCKER CONTAINERS ~~~~~ #
build: base gatk bwa fastqc trimmomatic sambamba bedtools variant-calling multiqc R msisensor delly2 annovar

pull:
	docker pull stevekm/ngs580-nf:base
	docker pull stevekm/ngs580-nf:bwa-0.7.17
	docker pull stevekm/ngs580-nf:fastqc-0.11.7
	docker pull stevekm/ngs580-nf:trimmomatic-0.36
	docker pull stevekm/ngs580-nf:sambamba-0.6.6
	docker pull stevekm/ngs580-nf:bedtools-2.26.0
	docker pull stevekm/ngs580-nf:variant-calling-0.0.2
	docker pull stevekm/ngs580-nf:multiqc-1.4
	docker pull stevekm/ngs580-nf:msisensor-0.2
	# docker pull broadinstitute/gatk:4.0.1.2
	# docker pull stevekm/ngs580-nf:R-3.2.3


# ~~~~~ CONVERT DOCKER CONTAINER TO SINGULARITY CONTAINER LOCALLY ~~~~~ #
REMOTE_CONTAINER_DIR:=/ifs/data/molecpathlab/containers/Singularity
USERNAME:=$(shell whoami)
SERVER:=phoenix.med.nyu.edu

docker2singularity: base bedtools bwa fastqc msisensor multiqc sambamba trimmomatic variant-calling R delly2
	for item in base bedtools-2.26.0 bwa-0.7.17 fastqc-0.11.7 msisensor-0.2 multiqc-1.4 sambamba-0.6.6 trimmomatic-0.36 variant-calling-0.0.2 R-3.3.2 delly2-0.7.7; \
	do ./docker2singularity.py "$$item"; \
	done

clean-singularity-imagefiles:
	find . -maxdepth 2 -type f -name "*.img" -delete

upload-imagefiles:
	rsync -vrthPlz -e ssh ./ $(USERNAME)@$(SERVER):$(REMOTE_CONTAINER_DIR) --include="*/" --include="*.img" --exclude="*" --prune-empty-dirs




# ~~~~~~ TEST CONTAINERS ~~~~~ #
test-base: base
	docker run --rm -ti stevekm/ngs580-nf:base bash

test-bwa: bwa
	docker run --rm -ti stevekm/ngs580-nf:bwa-0.7.17 bash

test-fastqc: fastqc
	docker run --rm -ti stevekm/ngs580-nf:fastqc-0.11.7 bash

test-gatk: gatk
	docker run --rm -ti broadinstitute/gatk bash

test-trimmomatic: trimmomatic
	docker run --rm -ti stevekm/ngs580-nf:trimmomatic-0.36 bash

test-sambamba: sambamba
	docker run --rm -ti stevekm/ngs580-nf:sambamba-0.6.6 bash

test-bedtools: bedtools
	docker run --rm -ti stevekm/ngs580-nf:bedtools-2.26.0 bash

test-variant-calling1: variant-calling1
	docker run --rm -ti stevekm/ngs580-nf:variant-calling-0.0.1 bash

test-variant-calling2: variant-calling2
	docker run --rm -ti stevekm/ngs580-nf:variant-calling-0.0.2 bash

test-variant-calling: test-variant-calling2

test-multiqc: multiqc
	docker run --rm -ti stevekm/ngs580-nf:multiqc-1.4 bash

test-msisensor: msisensor
	docker run --rm -ti stevekm/ngs580-nf:msisensor-0.2 bash

test-R: R
	docker run --rm -ti stevekm/ngs580-nf:R-3.3.2 bash

test-delly2: delly2
	docker run --rm -ti stevekm/ngs580-nf:delly2-0.7.7 bash

test-annovar: annovar
	docker run --rm -ti stevekm/ngs580-nf:annovar-150617 bash
