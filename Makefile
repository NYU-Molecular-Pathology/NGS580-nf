# Makefile to run the pipeline
SHELL:=/bin/bash
REFDIR:=/ifs/data/sequence/results/external/NYU/snuderllab/ref
EP:=
.PHONY: containers

# no default action
none:

# ~~~~~ SETUP PIPELINE ~~~~~ #
./nextflow:
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

bin/multiqc-venv/bin/activate:
	cd bin && \
	make -f multiqc.makefile setup

ref:
	[ -d "$(REFDIR)" ] && ln -fs $(REFDIR) ref || make -f ref.makefile ref

setup: install ref bin/multiqc-venv/bin/activate

containers:
	cd containers && make build

demo:
	git clone https://github.com/NYU-Molecular-Pathology/NGS580-demo-data.git

# ~~~~~ RUN PIPELINE ~~~~~ #
NGS580: setup ref
	./nextflow run main.nf  -with-dag flowchart-NGS580.dot $(EP) && \
	[ -f flowchart-NGS580.dot ] && dot flowchart-NGS580.dot -Tpng -o flowchart-NGS580.png

NGS580r: setup ref
	./nextflow run main.nf -resume -with-dag flowchart-NGS580.dot $(EP) && \
	[ -f flowchart-NGS580.dot ] && dot flowchart-NGS580.dot -Tpng -o flowchart-NGS580.png

NGS580l: install ref
	./nextflow run test.nf -profile local -with-dag flowchart-NGS580.dot $(EP) && \
	[ -f flowchart-NGS580.dot ] && dot flowchart-NGS580.dot -Tpng -o flowchart-NGS580.png

NGS580lr: install ref
	./nextflow run test.nf -profile local -resume -with-dag flowchart-NGS580.dot $(EP) && \
	[ -f flowchart-NGS580.dot ] && dot flowchart-NGS580.dot -Tpng -o flowchart-NGS580.png


# ~~~~~ CLEANUP ~~~~~ #
clean-traces:
	rm -f trace*.txt.*

clean-logs:
	rm -f .nextflow.log.*

clean-reports:
	rm -f *.html.*

clean-flowcharts:
	rm -f *.dot.*

clean-output:
	[ -d output ] && mv output oldoutput && rm -rf oldoutput &

clean-work:
	[ -d work ] && mv work oldwork && rm -rf oldwork &

clean: clean-logs clean-traces clean-reports clean-flowcharts

clean-all: clean clean-output clean-work
	[ -d .nextflow ] && mv .nextflow .nextflowold && rm -rf .nextflowold &
	rm -f .nextflow.log
	rm -f *.png
	rm -f trace*.txt*
	rm -f *.html*
