RMD_TEMPLATE:=main.Rmd
TIMESTAMP:=$(shell date +"%Y-%m-%d_%H-%M-%S")
none:

html:
	module load pandoc/1.13.1 && \
	compile_log="logs/compile_html.$(TIMESTAMP).log" && \
	{ \
	Rscript -e 'rmarkdown::render(input = "$(RMD_TEMPLATE)", params = list(input_dir = "/ifs/data/molecpathlab/development/NGS580-nf/output/analysis"))' ; \
	} 2>&1 | tee -a "$$compile_log" && \
	echo ">>> log saved to: $$compile_log"
	
	
