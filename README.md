# NGS580-nf
Target exome analysis for 580 gene panel

# Installation

- clone this repository

```
git clone https://github.com/NYU-Molecular-Pathology/NGS580-nf.git
cd NGS580-nf
```

- load Singularity / Docker containers

... [ coming soon ] ...

- download reference data

```
make ref
```

# Usage

- Create a samplesheet for your input `.fastq.gz` files:

```bash
./generate-samplesheets.py /path/to/fastq_dir
```

- [Optional] Manually create a [`samples.tumor.normal.csv` samplesheet](https://github.com/NYU-Molecular-Pathology/NGS580-nf/blob/master/example/samples.tumor.normal.csv) with the tumor-normal groupings for your samples, and update the original samplesheet with it by running the following script:

```bash
./update-samplesheets.py
```

- Run the pipeline:

```bash
make NGS580
```
__NOTE:__ This should be run in a `screen` session, or submitted to the HPC with the following script:

```bash
./submit.sh NGS580
```
