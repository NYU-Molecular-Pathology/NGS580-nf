#!/bin/bash

# submit the Nextflow pipeline as a job to run on the SGE HPC
# NOTE: requires that all nodes can qsub jobs
args="$@"
qsub_logdir="logs"
mkdir -p "${qsub_logdir}"

job_name="NGS580-nf"

qsub -wd $PWD -o :${qsub_logdir}/ -e :${qsub_logdir}/ -j y -N "$job_name" -b y bash -c "set -x; $args"
