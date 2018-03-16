#!/bin/bash

# environment configuration on NYU phoenix

module purge
module load local
module load sge/2011.11p1

module load singularity/2.4.2
module load fastqc/0.11.7
module load bwa/0.7.17
module load samtools/1.3
module load bedtools/2.26.0
module load java/1.8
module load jre/1.8
module load python/2.7.3

export PATH="/ifs/data/sequence/results/external/NYU/snuderllab/bin/trimmomatic/0.33:${PATH}" # trimmomatic.sh
export PATH="/ifs/data/sequence/results/external/NYU/snuderllab/bin/sambamba:${PATH}" # sambamba -> sambamba_v0.6.7
export PATH="/ifs/data/molecpathlab/bin/GenomeAnalysisTK-3.8-0:${PATH}" # gatk.sh
export PATH="/ifs/data/sequence/results/external/NYU/snuderllab/bin/lofreq/lofreq_star-2.1.2/bin:${PATH}" # lofreq
export PATH="/ifs/data/molecpathlab/bin/msisensor:${PATH}" # msisensor
# export PATH="/ifs/home/kellys04/software/delly/src:${PATH}" # delly
# export PATH="/ifs/home/kellys04/software/delly/src/bcftools:${PATH}" # bcftools
