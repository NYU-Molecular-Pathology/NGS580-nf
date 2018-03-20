VERSION:=150617
URL:=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.revision$(VERSION).tar.gz
ANNOVAR_DIR:=/ifs/data/molecpathlab/bin/annovar

none:

# symlink if dir exists, or download specified version
annovar: 
	[ -d "$(ANNOVAR_DIR)" ] && { \
	ln -fs "$(ANNOVAR_DIR)" annovar ; \
	} || { \
	wget "$(URL)" && \
	tar xvfz "annovar.revision$(VERSION).tar.gz" && \
	rm -f "annovar.revision$(VERSION).tar.gz" && \
	mv annovar "annovar-$(VERSION)" && \
	ln -fs "annovar-$(VERSION)" annovar ; \
	}

# symlink to scripts used in pipeline
setup: annovar
	[ -d "$(ANNOVAR_DIR)" ] && { \
	ln -fs "$(ANNOVAR_DIR)"/annotate_variation.pl ; \
	ln -fs "$(ANNOVAR_DIR)"/convert2annovar.pl ; \
	} || { \
	ln -fs "annovar-$(VERSION)"/annotate_variation.pl ; \
	ln -fs "annovar-$(VERSION)"/convert2annovar.pl ; \
	}

annotate_variation.pl convert2annovar.pl: setup

install: annotate_variation.pl convert2annovar.pl

clean-annovar:
	[ -d "annovar-$(VERSION)" ] && mv "annovar-$(VERSION)" oldannovar && rm -rf oldannovar &
	[ -L table_annovar.pl ] && rm -f table_annovar.pl || :
	[ -L coding_change.pl ] && rm -f coding_change.pl || :
	[ -L convert2annovar.pl ] && rm -f convert2annovar.pl || :
	[ -L retrieve_seq_from_fasta.pl ] && rm -f retrieve_seq_from_fasta.pl || :
	[ -L table_annovar.pl ] && rm -f table_annovar.pl || :
	[ -L variants_reduction.pl ] && rm -f variants_reduction.pl || :
	[ -L annotate_variation.pl ] && rm -f annotate_variation.pl || :
	[ -L annovar ] && unlink annovar
	rm -f annovar*.tar.gz
