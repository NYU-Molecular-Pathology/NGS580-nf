#!/bin/bash
# script to convert Docker container to Singularity image file

dir_tag="$1"
docker run -v /var/run/docker.sock:/var/run/docker.sock -v ${PWD}/${dir_tag}:/output --privileged -t --rm singularityware/docker2singularity stevekm/ngs580-nf:${dir_tag}
