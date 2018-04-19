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

- `Docker.makefile`: shortcuts to common actions for building Docker containers in this directory

- `Singularity.makefile`: shortcuts to common actions for building Singularity images in this container, and converting Docker images to Singularity container image files.

- `singularity.images.txt`: file containing names of Singularity image files that have already been built and transferred to the remote server for use

# Usage

## Build Docker Containers

You can build all Docker containers used in the pipeline with this command:

```
make build-all-Docker
```
- __NOTE:__ This takes up a lot of disk space. Make sure your Docker is set for a global image size of >100GB. It is also recommended to set Docker to use at least 8GB of memory when running. 

To build a selected container, run a command in the following format:

```
make -f Docker.makefile build VAR=<container_subdir>
```

where `container_subdir` is a subdirectory in this directory.

## Pull Docker Containers

Alternatively, you can pull these containers from Docker Hub with this command:

```
make pull-all-Docker
```

- __NOTE:__ Due to inconsistencies with Dockerhub, building is the most reliable way to ensure that you have the most recent version present on your system.

## Test a Docker Container

You can test out a Docker container by running the following command:

```
make -f Docker.makefile test VAR=<container_subdir>
```

This will start an interactive session inside the container, for you verify that the container works and the desired programs are installed.

## Converting Docker containers to Singularity containers

Once all Docker containers are present on your system, they can be converted to Singularity container image files for use on HPC systems.

To avoid creation & transfer of duplicate containers to the remote system, the `singularity.images.txt` file can be used to hold a list of known containers that already exist on the remote system.


The following command can be used to selectively convert only the missing containers from Docker format to Singularity format:

```
make convert-all-Docker
```
- __NOTE:__ to convert all containers to Singularity, simply clear the contents of the `singularity.images.txt` file, or run the Makefile command for the desired container in the format: `make -f Singularity.makefile convert-Docker VAR=<container_subdir>`.

Once built, all image files present in the current directory tree can be transferred to the remote HPC server using the following command:
```
make -f Singularity.makefile upload-imagefiles
```
- __NOTE__: see configuration details inside the `Singularity.makefile` file for determining the remote server address and directory path, which can be passed as Makefile command line arguments in the format `USERNAME=steve SERVER=remote.edu`, etc.

After uploading of images is finished, you can update the local copy of `singularity.images.txt` with the following command:

```
make -f Singularity.makefile update-remote-imagefile-list
```

## Test a Singularity image file

After generating a `.img` Singularity container image file, you can test it out by running the following command:

```
make -f Singularity.makefile test IMG=path/to/container.img
```

- __NOTE__: Vagrant is required for this. Tested on macOS only. See Vagrant installation notes in `VAGRANT.md`.

- __NOTE__: The relative path to the image file from the current directory should be passed to the `IMG` variable, not the absolute path.

This will load a copy of the official `singularityware` Vagrant image, bind in the path to your Singularity container file, and start an interactive Singularity shell session for you to use for verification that your programs are installed correctly inside the container.
