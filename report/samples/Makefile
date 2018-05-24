RMD_TEMPLATE:=main.Rmd
TIMESTAMP:=$(shell date +"%Y-%m-%d_%H-%M-%S")
FMT:=html_document
OUTPUTFILE:=report
INPUTDIR:=input
SAMPLEID:=
# all, html_document, pdf_document

compile:
	if [ "$$(module 2>/dev/null; echo $$? )" -eq "0" ]; then \
	module load pandoc/1.13.1 ; \
	module unload r ; \
	module load r/3.3.0 ; \
	fi ; \
	compile_log="logs/compile_html.$(TIMESTAMP).log" && \
	mkdir -p logs && \
	ln -fs main.Rmd "$(OUTPUTFILE).Rmd" &&  \
	{ \
	Rscript -e 'rmarkdown::render(input = "$(OUTPUTFILE).Rmd", params = list(input_dir = "$(INPUTDIR)", sampleID ="$(SAMPLEID)"), output_format = "$(FMT)", output_dir = ".")' ; \
	} 2>&1 | tee -a "$$compile_log" && \
	echo ">>> log saved to: $$compile_log"

pull:
	rsync -vrthP -e ssh kellys04@phoenix.med.nyu.edu:/ifs/data/molecpathlab/development/NGS580-nf/output/analysis output --exclude='*.bam' --exclude='*.fastq.gz' --exclude='*.sam'

clean:
	rm -f *.html
	rm -f *.Rdata
	rm -f *.pdf
	rm -f *.md
	rm -f *.tex
	[ -d main_files ] && rm -rf main_files
	
