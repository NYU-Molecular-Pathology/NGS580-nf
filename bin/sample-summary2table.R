#!/usr/bin/env Rscript

# this script will make a summary of .sample_summary files produced by GATK DepthOfCoverage

# ~~~~~~ FUNCTIONS ~~~~~~~ #
read.sample_summary <- function(file){
    # file <- input_file
    # read a GATK .csv formatted DepthOfCoverage sample_summary file into a dataframe
    df <- read.delim(file = file, header = TRUE, sep = ',', row.names = 1, check.names = FALSE)
    
    # add 'sample' column from rownames
    df[['Sample']] <- rownames(df)
    
    # reorder columns
    df <- df[c("Sample", colnames(df)[which(colnames(df) != 'Sample')])]
    
    # rename some columns
    names(df)[names(df) == 'mean'] <- 'MeanCoverage'
    names(df)[names(df) == 'granular_median'] <- 'MedianCoverage'
    names(df) <- gsub(pattern = '%_bases_above_', replacement = '', x = names(df))
    
    # remove some columns
    drops <- c("total","granular_third_quartile", "granular_first_quartile")
    df <- df[ , !(names(df) %in% drops)]
    
    # remove the 'Total' row
    df <- df[which(rownames(df) != 'Total'), ]
    rownames(df) <- c()
    return(df)
}



# ~~~~~~ RUN ~~~~~~ #
args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]

sample_summary_df <- read.sample_summary(input_file)

write.table(x = sample_summary_df, file = output_file, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# save.image("loaded.Rdata")
