# NGS580-nf

Target exome analysis for 580 gene panel (NGS580)

__NOTE:__ Details listed here may change during development

## Overview

This pipeline is designed to run targeted exome analysis on Illumina Next-Gen sequencing genomic data, in support of the NGS580 cancer diagnostic panel for NYU's Molecular Pathology Department.

This pipeline starts from paired-end fastq data (`.fastq.gz`), and is meant to accompany the output from the Illumina demultiplexing pipeline listed here: https://github.com/NYU-Molecular-Pathology/demux-nf.

The NGS580-nf analysis workflow includes read trimming, QC, alignment, variant calling, annotation, and reporting, along with many other steps.

## Contents

Some key pipeline components included in this repository:

- `bin`: directory of custom scripts used throughout the pipeline

- `containers`: directory of container recipes (Docker, Singularity) for use with the pipeline

- `example`: directory of example samplesheets, etc., to show the format used with this pipeline

- `targets`: directory of target region .bed files included with the pipeline for typical analyses

- `Makefile`: A Makefile with recipes for configuring, starting, and managing the pipeline. This is meant to be the main interface between the end-user and the pipeline. The Makefile should be reviewed as-needed to familiarize yourself with the methods and configurations that are meant to be used for running and managing the pipeline.

- `main.nf`: the main Nextflow pipeline script

- `nextflow.config`: configuration file for the main Nextflow pipeline script

- `.config.json`: a template for the required `config.json` file used in the pipeline, shows the default pipeline settings that are meant to be easily modified by the end-user and used within the pipeline for data processing.

- `annovar_db.nf`, `cnv-pool.nf`, `hapmap-pool.nf`, `ref.nf`: workflows for generating and downloading extra reference files used in the main pipeline.

- `ref`: default location for the storage of reference files (not used on NYU Big Purple HPC)

## Pipeline Items

Some key components that are created during setup, configuration, and execution of the pipeline:

- `samples.analysis.tsv`: the main samplesheet definig input items for the pipeline (described below)

- `config.json`: configuration file used for pipeline settings (see `.config.json` template for example)

- `output`: analysis output files published by the Nextflow pipeline

- `work`: Nextflow temporary directories for execution of pipeline tasks

- `trace.txt`, `nextflow.html`, `timeline.html`, `.nextflow.log`: Nextflow execution logs and reports

- `logs`: directory for pipeline execution logs

# Setup

This repository should first be cloned from GitHub:

```
git clone --recursive https://github.com/NYU-Molecular-Pathology/NGS580-nf.git
cd NGS580-nf
```

- Once a copy of the repo is made, it can be used to "deploy" new copies of the workflow in a pre-configured state

## Reference Data

Nextflow pipelines have been included for downloading required reference data, including ANNOVAR reference databases. You can run them with the following command:

```
make setup
```

### HapMap Pool .bam

A negative control HapMap pool .bam file can be prepared using the following command:

```
make hapmap-pool
```

- Requires `samples.hapmap.tsv` file specifying the .bam files to be combined (example included at `example/samples.hapmap.tsv`).

This file is typically built from multiple HapMap samples previously aligned by this pipeline. For demonstration purposes, you can provide any .bam and .bai files. These should be listed under the `HapMapBam` and `HapMapBai` keys of `config.json`.

### CNV Pool

A control normal sample .cnn file for CNV calling can be prepared using the following command:

```
make cnv-pool
```

- Requires `samples.cnv.tsv` file specifying the .bam files to be used (example included at `example/samples.cnv.tsv`)

This file is typically built from .bam files of specially chosen normal tissue sequencing samples previously aligned by this pipeline. For demonstration purposes, you can create the .cnn file from any desired .bam file. Note that the targets .bed file used to create the .cnn file must match the targets used in the rest of the pipeline. The .cnn file used should be set under the `CNVPool` key in `config.json`.

## Containers

The `containers` directory contains instructions and recipes for building the Docker and Singularity containers used in the pipeline.

Docker is typically used for local container development, while Singularity containers are used on the NYU Big Purple HPC cluster. The current pipeline configuration for Big Purple uses `.simg` files stored in a common location on the file system.

# Usage

The pipeline is designed to start from demultiplexed paired end `.fastq.gz` files, with sample ID, tumor ID, and matched normal ID associations defined for each set of R1 and R2 .fastq file using a file `samples.analysis.tsv` (example included at `example/samples.analysis.tsv`).

## Deployment

The easiset way to use the pipeline is to "deploy" a new instance of it based on output from the demultiplexing pipeline [`demux-nf`](https://github.com/NYU-Molecular-Pathology/demux-nf). This will automatically propagate configurations and information from the demultiplexing output.

The pipeline can also deploy a new, pre-configured copy of itself using the included `deploy` recipe:

```
make deploy PRODDIR=/path/to/NGS580_analyses RUNID=Name_for_analysis FASTQDIR=/path/to/fastq_files
```

- An optional argument `DEMUX_SAMPLESHEET` can be used to provide a specially formatted demultiplexing samplesheet to be used for extracting extra sample information (example included at `example/demux-SampleSheet.csv`; note the extra columns labeling tumor-normal pair IDs, used later).

## Create Config

A file `config.json` is required to hold settings for the pipeline. It should be created using the built-in methods:

```
make config RUNID=my_run_ID FASTQDIR=/path/to/fastqs
```

or

```
make config RUNID=my_run_ID FASTQDIRS='/path/to/fastqs1 /path/to/fastqs2'
```

Once created, the `config.json` file can be updated manually as needed. The template and default values can be viewed in the included `.config.json` file.

- `config.json` should be generated automatically if you used `make deploy`

## Create Samplesheet

A samplesheet file `samples.analysis.tsv` is required in order to define the input samples and their associated .fastq files (example included at `example/samples.analysis.tsv`). Create a samplesheet, based on the config file, using the built-in methods:

```
make samplesheet
```

- Note that this uses the values previously saved in `config.json` to create the samplesheet

### Sample Pairs (Optional)

The NGS580-nf pipeline has special processing for tumor-normal pairs. These pairs should be defined in the `samples.analysis.tsv` file, by listing the matched Normal sample for each applicable sample.

In order to update `samples.analysis.tsv` automatically with these sample pairs, an extra samplesheet can be provided with the tumor-normal pairs.

Create a `samples.tumor.normal.csv` samplesheet (example included at `example/samples.tumor.normal.csv`) with the tumor-normal groupings for your samples, and update the original samplesheet with it by running the following script:

```
python update-samplesheets.py --tumor-normal-sheet samples.tumor.normal.csv
```

If a demultiplexing samplesheet with extra tumor-normal pairs information was supplied (see example: `example/demux-SampleSheet.csv`), then it can be used to update the samplesheet with pairs information with the following recipe:

```
make pairs PAIRS_SHEET=demux-SampleSheet.csv PAIRS_MODE=demux
```

## Run

The pipeline includes an auto-run functionality that attempts to determine the best configuration to use for NYU phoenix and Big Purple HPC clusters:

```
make run
```

This will run the pipeline in the current session.

In order to run the pipeline in the background as a job on NYU's Big Purple HPC, you should instead use the `submit` recipe:

```
make submit SUBQ=fn_medium
```

Where `SUBQ` is the name of the SLURM queue you wish to use.

Refer to the `Makefile` for more run options.

Due to the scale of the pipeline, a "local" run option is not currently configured, but can be set up easily based on the details shown in the Makefile and `nextflow.config`.

### Extra Parameters

You can supply extra parameters for Nextflow by using the `EP` variable included in the Makefile, like this:

```
make run EP='--runID 180320_NB501073_0037_AH55F3BGX5
```

### Demo

A demo dataset can be loaded using the following command:

```
make demo
```

This will:

- checkout a demo dataset

- create a `samples.analysis.tsv` samplesheet for the analysis

You can then proceed to run the analysis with the commands described above.

## More Functionality

Extra functions included in the Makefile for pipeline management include:

### `make clean`

Removes all Nextflow output except for the most recent run. Use `make clean-all` to remove all pipeline outputs.

### `make record PRE=some_prefix_`

"Records" copies of the most recent pipeline run's output logs, configuration, Nextflow reports, etc.. Useful for recording analyses that failed or had errors in order to debug. Include the optional argument `TASK` to specify a Nextflow `work` directory to include in the records (example: `make record PRE=error_something_broke_ TASK=e9/d9ff34`).

### `make kill`

Attempts to cleanly shut down a pipeline running on a remote host e.g. inside a SLURM HPC compute job. Note that you can also use `scancel` to halt the parent Nextflow pipeline job as well.

### `make fix-permissions`, `make fix-group`

Attempts to fix usergroup and permissions issues that may arise on shared systems with multiple users. Be sure to use the extra argument `USERGROUP=somegroup` to specify the usergroup to update to.

### `make finalize-work-rm`

Examines the `trace.txt` output from the most recent completed pipeline run in order to determine while subdirectories in the Nextflow `work` dir are no longer needed, and then deletes them. Can delete multiple subdirs in parallel when run with `make finalize-work-rm -j 20` e.g. specifying to delete 20 at a time, etc.

# Software

Developed under Centos 6, RHEL 7, macOS 10.12

- bash

- GNU `make`, standard GNU tools

- Python 2/3

- Java 8+ for Nextflow

- Docker/Singularity as needed for containers
