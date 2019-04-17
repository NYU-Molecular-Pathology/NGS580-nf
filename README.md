# NGS580-nf
Target exome analysis for 580 gene panel

__NOTE:__ Details listed here may change during development

# Installation

- clone this repository

```
git clone https://github.com/NYU-Molecular-Pathology/NGS580-nf.git
cd NGS580-nf
```

## Docker/Singularity

See instructions in the `containers` directory for building Docker and Singularity images. At this time, we are only using Singularity images on NYU Big Purple HPC, so those should be considered the 'production' containers while the Docker containers are used more for development and testing purposes.


## Reference Data

Nextflow pipelines have been included for downloading required reference data, including ANNOVAR reference databases. You can run them with the following command:

```
make setup
```

Note that this will require a version of ANNOVAR (included Docker container is used by default), along with Nextflow, which requires Java 8 and will be automatically installed in the current directory.

# Usage

The pipeline is designed to start from demultiplexed `.fastq.gz` files, produced by an Illumina sequencer.

## Create Config

A config file (`config.json`) should first be created using the built-in methods:

```
make config RUNID=my_run_ID FASTQDIR=/path/to/fastqs
```

or

```
make config RUNID=my_run_ID FASTQDIRS='/path/to/fastqs1 /path/to/fastqs2'
```

Once created, the `config.json` file can be updated manually as needed.

## Create Samplesheet

Create a samplesheet, based on the config file, using the built-in methods:

```
make samplesheet
```

- [Optional] Manually create a [`samples.tumor.normal.csv` samplesheet](https://github.com/NYU-Molecular-Pathology/NGS580-nf/blob/master/example/samples.tumor.normal.csv) with the tumor-normal groupings for your samples, and update the original samplesheet with it by running the following script:

```
./update-samplesheets.py
```

## Run the pipeline:

If you are on NYU's phoenix or Big Purple HPC clusters, you can use the auto-run functionality based on the config file:

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

### Extra Parameters

You can supply extra parameters for Nextflow by using the `EP` variable included in the Makefile, like this:

```
make run EP='--runID 180320_NB501073_0037_AH55F3BGX5 
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
