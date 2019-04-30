# HapMap Pool

Synthetic 'normal' human DNA sample, representing typical background SNP's. Many sequencing runs include a HapMap sample as a negative control to simulate a sample that should not have any clinically relevant variants.

A pool of HapMap samples can be used for quality control, by performing variant calling against the HapMap sample in a given sequencing run against a pooled HapMap sample, where the current run's HapMap is treated as the 'tumor' sample and the pool is treated as the 'normal' sample. Any variants called this way would represent potential contaminants, since the HapMap should not have any new variants in it.

## Creating a HapMap Pool

First, set up the pipeline with .fastq files as normal, for each HapMap sample you wish to include in the pool. Run the full standard pipeline with these samples.

Next, create the file `samples.hapmap.tsv` in the main pipeline directory. It should be a tab-separated file with headers formatted like this:

```
Sample	Bam	Bai
Run1-HapMap output/alignments/recalibrated/Run1-HapMap.bam  output/alignments/recalibrated/Run1-HapMap.bam.bai
Run2-HapMap output/alignments/recalibrated/Run2-HapMap.bam  output/alignments/recalibrated/Run2-HapMap.bam.bai
```

Then, run the specialized `hapmap-pool.nf` pipeline:

```
make hapmap-pool
```

or

```
nextflow run hapmap-pool.nf -profile hapmap_pool
```

Note that the `hapmap_pool` profile in the `nextflow.config` is currently configured to run on NYUMC's Big Purple HPC system, and may need to be modified to run on your system.

Running with 40 CPU cores and 6 HapMap samples (19GB total) takes about 1.5 hours.

The output will be `output-hapmap-pool/HapMap-pool.bam` and `output-hapmap-pool/HapMap-pool.bam.bai`.

## Resources

Details about HapMap sample can be found here: https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Product=CC&Ref=GM12878
