VAR=

none:

.PHONY: base $(VAR)

check-Docker:
	docker --version > /dev/null 2>&1 || { echo "ERROR: 'docker' not found" && exit 1 ; }

# ~~~~~ BUILD DOCKER CONTAINERS ~~~~~ #
base: check-Docker
	cd base && docker build -t stevekm/ngs580-nf:base .

# ~~~~~~~ SETUP DOCKER CONTAINERS ~~~~~ #
build: base
	cd $(VAR) && \
	docker build -t stevekm/ngs580-nf:$(VAR) .

pull:
	docker pull stevekm/ngs580-nf:$(VAR)

# ~~~~~~ TEST CONTAINERS ~~~~~ #
test-base: base
	docker run --rm -ti stevekm/ngs580-nf:base bash

test: build
	docker run --rm -ti stevekm/ngs580-nf:$(VAR) bash
