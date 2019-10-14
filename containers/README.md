# NGS580-nf Containers

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/744)

Docker and Singularity containers for use with the pipeline.

Pre-built copies of these containers are hosted here:

- [Docker Hub](https://hub.docker.com/r/stevekm/ngs580-nf/)

- [Singularity Hub](https://www.singularity-hub.org/collections/744)

- NOTE: pre-built copies on Docker Hub and Singularity Hub may not be up to date, for best results build directly from the included Dockerfiles and Singularity recipe files.

Instructions here are meant to be used on a local desktop system (e.g. macOS). If you are running Linux, some of these commands will not be necessary.

# Contents

Items in this directory:

- `Makefile`: shortcuts to common actions for using and configuring containers

- `Vagrantfile`: configuration for Vagrant virtual machine to run Singularity for building Singularity images (instructions for installing Vagrant in `VAGRANT.md`)

- [other directories]: Dockerfile and Singularity files for other tools used in the pipeline

# Usage

## Docker

Build a Docker container in this directory with a command like this:

```
make docker-build VAR=<dirname>
```

example:

```
make docker-build VAR=htslib-1.7
```

Test a Docker image with the command:

```
make docker-test VAR=<dirname>
```

## Singularity

There are two methods to build Singularity images on macOS; running Singularity inside Docker, or running Singularity inside Vagrant virtual machine.

When using the include Vagrant methods, the Vagrant virtual machine with Singularity inside will be built automatically.

Build a Singularity container in this directory with Vagrant:

```
make singularity-build VAR=<dirname>
```

Test the Singularity container with Vagrant:

```
make singularity-test VAR=<dirname>
```

To build with Docker instead (required for some containers, especially those that bootstrap from other Docker containers), first build the Singularity Docker image:

```
make docker-build VAR=Singularity-2.4
```

This container will then be used to build Singularity containers inside Docker with commands like:

```
make singularity-build-docker VAR=<dirname>
```

The final step for using Singularity containers with the pipeline is to prepare to copy the output `.simg` file over to the HPC cluster, using a command like this:

```
make singularity-sync VAR=<dirname>
```

This will print to console an example `rsync` command to be used to copy the image file over. Review the command, and remove the `--dry-run` arg to perform the copy operation after you are satisfied that the `rsync` command is correct. 

# Software

Tested with

- Docker version 17.12.0-ce, build c97c6d6

- Vagrant 2.0.1

- Singularity 2.4
