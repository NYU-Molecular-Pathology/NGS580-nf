# NGS580-nf
Target exome analysis for 580 gene panel

__NOTE:__ Details listed here may change during development

# Installation 

- clone this repository

```
git clone https://github.com/NYU-Molecular-Pathology/NGS580-nf.git
cd NGS580-nf
```

# Setup 

## Docker

If you are using the included Docker containers (recommended), run the following commands to build the included containers:

```
cd containers
make build-all-Docker
```

## Conda

Conda environment recipes (equivalent to the included Docker containers) have also been included and can be created with the following commands:

```
cd conda
make 
```

## Reference Data

Nextflow pipelines have been included for downloading required reference data, including ANNOVAR reference databases. You can run them with the following command:

```
make setup
```

Note that this will require a version of ANNOVAR (included Docker container is used by default), along with Nextflow, which requires Java 8 and will be automatically installed in the current directory.

# Usage

The pipeline is designed to start from demultiplexed .fastq.gz files, produced by an Illumina sequencer, which have been split by flowcells lane (default `bcl2fastq` output format).

## Create Samplesheet

A samplesheet is needed to gather data about fastq.gz input files and sample IDs. It can be generated with these steps:

- Create a `samples.analysis.tsv` samplesheet (example [here](https://github.com/NYU-Molecular-Pathology/NGS580-nf/blob/master/example/samples.analysis.tsv)) for your input `.fastq.gz` files:

```
./generate-samplesheets.py /path/to/fastq_dir
```

- [Optional] Manually create a [`samples.tumor.normal.csv` samplesheet](https://github.com/NYU-Molecular-Pathology/NGS580-nf/blob/master/example/samples.tumor.normal.csv) with the tumor-normal groupings for your samples, and update the original samplesheet with it by running the following script:

```
./update-samplesheets.py
```

## Run the pipeline:

To run the Nextflow pipeline in the current session:

- on NYU phoenix HPC using 'module', submits jobs to SGE scheduler
    
```
make run-phoenix
```

- on MCIT Power8 server using Miniconda, runs jobs on host system

```
make run-power

```

- locally using Docker, runs jobs on host system

```
make run-local
```

To submit the parent Nextflow process as a job on the HPC system:

- on NYU phoenix HPC

```
make submit-phoenix
```

# Demo

A demo dataset can be loaded using the following command:

```
make demo
```

This will:

- checkout a demo dataset

- create a `samples.analysis.tsv` samplesheet for the analysis

You can then proceed to run the analysis with the commands described above.
