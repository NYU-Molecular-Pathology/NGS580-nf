# setup for MultiQC virtual environment in Python 2.7 on phoenix HPC
none:

multiqc-venv/bin/activate:
	module unload python && \
	module load python/2.7.3 && \
	export PYTHONPATH= && \
	virtualenv multiqc-venv --no-site-packages

multiqc-activate: multiqc-venv/bin/activate
	ln -fs multiqc-venv/bin/activate multiqc-activate

install: multiqc-venv/bin/activate multiqc-activate
	export PYTHONPATH= && \
	source multiqc-venv/bin/activate && \
	pip install -r multiqc.requirements.txt

setup: install multiqc-activate

clean:
	[ -d multiqc-venv ] && mv multiqc-venv multiqc-venvold && rm -rf multiqc-venvold &
	unlink multiqc-activate

test: multiqc-venv/bin/activate
	export PYTHONPATH= && \
	source multiqc-venv/bin/activate && \
	multiqc --version