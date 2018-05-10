#!/usr/bin/env Rscript

# this script will make a consistent table of .sample_interval_summary files produced by GATK DepthOfCoverage for usage elsewhere in the pipeline

# ~~~~~~ FUNCTIONS ~~~~~~~ #
chrom_regions2df <- function(regions){
    # split the regions into chrom coordinates for BED files
    # regions <- c("chr1:236998847-236998987", "chr1:237001714-237001899")
    # regions <- df[['Target']]
    regions_df <- as.data.frame(do.call(rbind, strsplit(regions, ':')))
    regions_df <- cbind(regions_df[1],
                        as.data.frame(do.call(rbind, strsplit(as.character(regions_df$V2), '-'))))
    colnames(regions_df) <- c("Chrom", "Start", "Stop")
    return(regions_df)
}

read.sample_interval_summary <- function(file){
    # file <- input_file
    df <- read.delim(file = file, header = TRUE, sep = ',', check.names = FALSE, stringsAsFactors = FALSE)
    
    # remove some columns
    drops <- c("total_coverage","average_coverage")
    df <- df[ , !(names(df) %in% drops)]
    df <- df[, names(df)[!grepl(pattern = '_granular_Q1', x = names(df))]]
    df <- df[, names(df)[!grepl(pattern = '_granular_Q3', x = names(df))]]
    
    # rename some columns
    names(df)[grepl(pattern = '_total_cvg', x = names(df))] <- 'TotalCoverage'
    names(df)[grepl(pattern = '_mean_cvg', x = names(df))] <- 'MeanCoverage'
    names(df)[grepl(pattern = '_granular_median', x = names(df))] <- 'MedianCoverage'
    names(df) <- gsub(pattern = '^.*%_above_', replacement = '', x = names(df), perl = TRUE)
    
    # add regions columns
    df <- cbind(chrom_regions2df(df[['Target']]), df)
    
    return(df)
}

# ~~~~~~ RUN ~~~~~~ #
args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]

sample_interval_summary_df <- read.sample_interval_summary(input_file)

write.table(x = sample_interval_summary_df, file = output_file, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# save.image("loaded.Rdata")
