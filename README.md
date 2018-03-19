# NGS580-nf
Target exome analysis for 580 gene panel

# Installation & Setup

- clone this repository

```
git clone https://github.com/NYU-Molecular-Pathology/NGS580-nf.git
cd NGS580-nf
```

## Setup Dependencies

### With Docker

If you are using Docker containers (recommended), run the following command to build the included containers:

```
make setup
```

This will:

- build Singularity / Docker containers

- download genomic reference data and ANNOVAR databases if the pre-configured data directories are not available; if they exist, they will be symlinked instead (see Makefile for details)

### NYU phoenix HPC

If running on the NYU phoenix HPC cluster, use the following command:

```
make setup-p
```

This will:

- set up reference data as per the `setup` command

- build a Python virtual environment for `multiqc` in the `bin` directory


# Usage

To run the pipeline, follow these steps:

- Create a `samples.analysis.tsv` samplesheet (example [here](https://github.com/NYU-Molecular-Pathology/NGS580-nf/blob/master/example/samples.analysis.tsv)) for your input `.fastq.gz` files:

```
./generate-samplesheets.py /path/to/fastq_dir
```

- [Optional] Manually create a [`samples.tumor.normal.csv` samplesheet](https://github.com/NYU-Molecular-Pathology/NGS580-nf/blob/master/example/samples.tumor.normal.csv) with the tumor-normal groupings for your samples, and update the original samplesheet with it by running the following script:

```
./update-samplesheets.py
```

- Run the pipeline:

```
make run
```
__NOTE:__ This should be run in a `screen` session, or submitted to the HPC with the following SGE submission script:

```
./submit.sh run
```

## Demo

A demo dataset can be loaded using the following command:

```
make demo
```

This will:

- checkout a demo dataset

- create a `samples.analysis.tsv` samplesheet for the analysis

You can then proceed to run the analysis with the commands described above.
