SHELL:=/bin/bash
VAR=
REMOTE_CONTAINER_DIR:=/ifs/data/molecpathlab/containers/Singularity
USERNAME:=$(shell whoami)
SERVER:=phoenix.med.nyu.edu
IMG=

none:

.PHONY: base $(VAR)

check-docker:
	echo ">>> Making sure Docker is present..."
	docker --version > /dev/null 2>&1 || { echo "ERROR: 'docker' not found" && exit 1 ; }

check-docker-image: check-docker
	echo ">>> Making sure Docker image exists..."
	[ "$$(docker images -q stevekm/ngs580-nf:$(VAR))" != "" ] && exit 0 || exit 1

convert-Docker: check-docker-image
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
	rm -f Vagrantfile*
	[ -d .vagrant ] && rm -rf .vagrant || :

test: check-vagrant clean-vagrant
	[ -z "$(IMG)" ] && { echo ">>> ERROR: No image file passed" ; exit 1 ; } || :
	[ ! -f "$(IMG)" ] && { echo ">>> ERROR: Invalid image file" ; exit 1 ; } || :
	echo ">>> Testing Singularity image file: $(IMG)"
	image_dir="$$(dirname $(IMG))" && \
	timestamp="$$(date +"%s")" && \
	vagrant init singularityware/singularity-2.4 && \
	tmpfile="Vagrantfile.$${timestamp}" && \
	echo ">>> Moving Vagrantfile to $${tmpfile}" && \
	/bin/mv Vagrantfile "$${tmpfile}" && \
	echo ">>> Writing new Vagrantfile" && \
	cat "$${tmpfile}" | sed -e "s|  # config.vm.synced_folder \"../data\", \"/vagrant_data\"|  config.vm.synced_folder \"$${image_dir}\", \"/$${image_dir}\"|" > Vagrantfile && \
	echo ">>> Starting Vagrant to test image file..." && \
	vagrant up && \
	vagrant ssh -c "singularity shell /$(IMG)"
