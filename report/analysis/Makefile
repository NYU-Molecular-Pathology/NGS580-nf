RMD_TEMPLATE:=main.Rmd
TIMESTAMP:=$(shell date +"%Y-%m-%d_%H-%M-%S")
none:

html:
	[ "$$(module 2>/dev/null; echo $$? )" -eq "0" ] && module load pandoc/1.13.1 || : ; \
	compile_log="logs/compile_html.$(TIMESTAMP).log" && \
	mkdir -p logs && \
	{ \
	Rscript -e 'rmarkdown::render(input = "$(RMD_TEMPLATE)", params = list(input_dir = "input"))' ; \
	} 2>&1 | tee -a "$$compile_log" && \
	echo ">>> log saved to: $$compile_log"
	
	
pull:
	rsync -vrthP -e ssh kellys04@phoenix.med.nyu.edu:/ifs/data/molecpathlab/development/NGS580-nf/output/analysis output --exclude='*.bam' --exclude='*.fastq.gz' --exclude='*.sam'
