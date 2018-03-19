# NGS580-nf Containers

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/744)

Docker and Singularity containers for use with the pipeline. 

Pre-built copies of these containers are hosted here:

- [Docker Hub](https://hub.docker.com/r/stevekm/ngs580-nf/)

- [Singularity Hub](https://www.singularity-hub.org/collections/744)

# Contents

Items in this directory:

- `base`: directory with Dockerfile and Singularity file for the base layer used for most images in this repository. Contains custom mount points and configurations. 

- `variant-calling-x.x.x`: directories with Dockerfile and Singularity file for containers with variant calling tools, built from Broad Institute's official GATK Docker containers.

- [other directories]: Dockerfile and Singularity files for other tools used in the pipeline

- `Makefile`: shortcuts to common actions for using and configuring containers

- `docker2singularity.sh`: simple script for converting a Docker container in this directory to a Singularity container

- `docker2singularity.py`: script that converts Docker containers to Singularity containers if the directory is not listed as already having an image file in the `singularity.images.txt` file

- `singularity.images.txt`: file containing names of Singularity image files that have already been built and transfered to the remote server for use

# Usage

## Build Docker Containers

You can build all Docker containers used in the pipeline with this command:

```
make build
```

## Pull Docker Containers

Alternatively, you can pull these containers from Docker Hub with this command:

```
make pull
```

## Converting Docker containers to Singularity containers

The Docker containers in this directory can be converted to Singularity containers for use on HPC systems with the following command:

```
make docker2singularity
```

Once built, these new image files can be transfered to the remote HPC server using the following command:
```
make upload-imagefiles
```
- NOTE: see configuration details inside the Makefile for determining the remote server address and directory path
