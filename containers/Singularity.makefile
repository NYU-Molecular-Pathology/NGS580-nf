SHELL:=/bin/bash
VAR=
REMOTE_CONTAINER_DIR:=/ifs/data/molecpathlab/containers/Singularity
USERNAME:=$(shell whoami)
SERVER:=phoenix.med.nyu.edu
IMG=

none:

.PHONY: base $(VAR)

check-Docker:
	echo ">>> Making sure Docker is present..."
	docker --version > /dev/null 2>&1 || { echo "ERROR: 'docker' not found" && exit 1 ; }

check-Docker-image: check-Docker
	echo ">>> Making sure Docker image exists..."
	[ "$$(docker images -q stevekm/ngs580-nf:$(VAR))" != "" ] && exit 0 || exit 1

convert-Docker: check-Docker-image
	echo ">>> Checking if an entry already exists in singularity.images.txt..."
	[ "$$(grep -q '$(VAR)' singularity.images.txt; echo $$?)" -eq 0 ] && \
	echo ">>> A Singularity image file is already known for this image, skipping..." || { \
	echo ">>> Converting Docker image $(VAR) to Singularity image file" ; \
	docker run -v /var/run/docker.sock:/var/run/docker.sock -v $${PWD}/$(VAR):/output --privileged -t --rm singularityware/docker2singularity stevekm/ngs580-nf:$(VAR) ; \
	}

clean-imagefiles:
	find . -maxdepth 2 -type f -name "*.img" -delete

upload-imagefiles:
	rsync -vrthPlz -e ssh ./ $(USERNAME)@$(SERVER):$(REMOTE_CONTAINER_DIR) --include="*/" --include="*.img" --exclude="*" --prune-empty-dirs

check-remote-imagefiles:
	ssh $(USERNAME)@$(SERVER) 'cd $(REMOTE_CONTAINER_DIR) && ls -l *'

update-remote-imagefile-list:
	timestamp="$$(date +"%s")" && \
	tmpfile="tmp.$${timestamp}" && \
	echo ">>> Temp file list will be written to: $${tmpfile}" && \
	ssh $(USERNAME)@$(SERVER) 'cd $(REMOTE_CONTAINER_DIR) && find . -type f -name "*.img" | sort | cut -sd / -f 2-' > "$${tmpfile}" && \
	old_file="singularity.images.old.$${timestamp}.txt" && \
	echo ">>> singularity.images.txt will be moved to $${old_file}" && \
	/bin/mv singularity.images.txt "$${old_file}" && \
	echo ">>> moving $${tmpfile} to singularity.images.txt" && \
	/bin/mv "$${tmpfile}" "singularity.images.txt"

check-vagrant:
	echo ">>> Making sure Vagrant is present..."
	vagrant -v >/dev/null 2>&1 || { echo "ERROR: 'docker' not found" && exit 1 ; }

clean-vagrant:
	[ -d .vagrant ] && rm -rf .vagrant || :

#  test with Vagrant
test: check-vagrant
	[ -z "$(IMG)" ] && { echo ">>> ERROR: No image file passed (e.g. IMG=path/to/container.img)" ; exit 1 ; } || :
	[ ! -f "$(IMG)" ] && { echo ">>> ERROR: Invalid image file (IMG=$(IMG))" ; exit 1 ; } || :
	echo ">>> Testing Singularity image file: $(IMG)"
	echo ">>> Starting Vagrant..." && \
	vagrant up test && \
	vagrant ssh test -c "singularity shell /vagrant/$(IMG)"

# build with Vagrant
build: check-vagrant
	[ -z "$(VAR)" ] && { echo ">>> ERROR: No directory variable file passed" ; exit 1 ; } || :
	[ ! -d "$(VAR)" ] && { echo ">>> ERROR: Invalid directory variable" ; exit 1 ; } || :
	[ ! -f "$(VAR)/Singularity.$(VAR)" ] && { echo ">>> ERROR: Singularity file '$(VAR)/Singularity.$(VAR)' does not exist" ; exit 1 ; } || :
	echo ">>> Setting up to build Singularity image in directory: $(VAR)"
	image_file="stevekm_ngs580-nf_$(VAR).img" && \
	image_path="$(VAR)/$${image_file}" && \
	[ -f "$${image_path}" ] && { echo ">>> Removing previous image file: $${image_path}" ; rm -f "$${image_path}" ; } ; \
	echo ">>> Output file will be: $(VAR)/$${image_file}" && \
	vagrant up build && \
	vagrant ssh build -c "cd /vagrant/$(VAR) && sudo singularity build $${image_file} Singularity.$(VAR)" && \
	echo ">>> Output file: $(VAR)/$${image_file}"
