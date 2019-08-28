#!/bin/bash -x
#SBATCH -D /gpfs/data/molecpathlab/development/NGS580-development-runs/run-2
#SBATCH -o /gpfs/data/molecpathlab/development/NGS580-development-runs/run-2/logs/slurm-%j.log.1567000969.out
#SBATCH -J NGS580-run-2
#SBATCH -p intellispace
#SBATCH --time=5-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -c 8
#SBATCH --mem 48G
#SBATCH --export=HOSTNAME
touch .nextflow.submitted
get_pid(){ head -1 .nextflow.pid; }
rm_submit(){ echo ">>> trap: rm_submit" ; [ -e .nextflow.submitted ] && rm -f .nextflow.submitted || : ; }
wait_pid(){ local pid=$1 ; while kill -0 $pid; do echo waiting for process $pid to end ; sleep 1 ; done ; }
nxf_kill(){ rm_submit ; echo ">>> trap: nxf_kill" && pid=$(get_pid) && kill $pid && wait_pid $pid ; }
trap nxf_kill HUP
trap nxf_kill INT
trap nxf_kill EXIT
make submit-bigpurple-run TIMESTAMP=1567000969 
