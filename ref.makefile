# Setup the reference files for the pipeline
SHELL:=/bin/bash
none:

_ref:
	wget -r --no-parent -e robots=off -nH --cut-dirs=4 https://genome.med.nyu.edu/results/external/NYU/snuderllab/ref/
	find ref/ -type f -name "index.html*" -delete

ref:
	wget https://genome.med.nyu.edu/results/external/NYU/snuderllab/ref.tar.gz && \
	tar -vxzf ref.tar.gz

HG19_GENOME_FA_MD5:=c1ddcc5db31b657d167bea6d9ff354f9
ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa: ref
	wget -r --no-parent -e robots=off -nH --cut-dirs=4 https://genome.med.nyu.edu/results/external/NYU/snuderllab/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/ && \
	bin/compare_md5.sh ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa "$(HG19_GENOME_FA_MD5)"

clean-ref:
	[ -d ref ] && mv ref oldref && rm -rf oldref &
