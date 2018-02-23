SHELL:=/bin/bash
none:

# ~~~~~ NEXTFLOW PIPELINE ~~~~~ #
# NextFlow setup & run commands
./nextflow:
	curl -fsSL get.nextflow.io | bash

install: ./nextflow

setup: 
	./generate-samplesheets.py example-data/ && ./update-samplesheets.py


NGS580: install
	./nextflow run main.nf  -with-dag flowchart-NGS580.dot && \
	[ -f flowchart-NGS580.dot ] && dot flowchart-NGS580.dot -Tpng -o flowchart-NGS580.png

NGS580r: install
	./nextflow run main.nf -resume -with-dag flowchart-NGS580.dot && \
	[ -f flowchart-NGS580.dot ] && dot flowchart-NGS580.dot -Tpng -o flowchart-NGS580.png

# entire pipeline
run: NGS580


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