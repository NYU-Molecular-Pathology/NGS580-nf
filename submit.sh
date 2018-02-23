#!/bin/bash

# submit the Nextflow pipeline as a job to run on the SGE HPC
# NOTE: requires that all nodes can qsub jobs

target="${1:-'xxxxxxxxxxxxxxxxxxx'}"
[ "$xxxxxxxxxxxxxxxxxxx" == 'xxxxxxxxxxxxxxxxxxx' ] && printf "please specify a target to make" && exit 1 || :

qsub_logdir="qsub-logs"
mkdir -p "${qsub_logdir}"

job_name="NGS580-nf"

qsub -wd $PWD -o :${qsub_logdir}/ -e :${qsub_logdir}/ -j y -N "$job_name" <<E0F
    set -x
    make "${target}"
E0F
