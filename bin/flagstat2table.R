#!/usr/bin/env Rscript

# this script will read in a Samtools/Sambamba flagstat file and output a table of values
# udpate main parsing function to return new columns as needed
# https://www.biostars.org/p/12475/

# ~~~~~~ FUNCTIONS ~~~~~~~ #
read.flagstat <- function(file){
    # reads a samtools flagstat text file to parse values
    # file <- input_file
    file_lines <- readLines(file)
    
    # total number of alignments
    mapped_all_line <- grep(pattern = 'mapped (', x = file_lines, value = TRUE, fixed = TRUE)[1]
    mapped_all <- as.numeric(unlist(strsplit(x = mapped_all_line, split = ' '))[1])
    
    # the read is paired in sequencing, no matter whether it is mapped in a pair
    total_paired_line <- grep(pattern = ' paired in sequencing', x = file_lines, value = TRUE, fixed = TRUE)[1]
    total_paired <- as.numeric(unlist(strsplit(x = total_paired_line, split = ' '))[1])
    
    # same as the 'properly paired' line
    # mapped_secondary_line <- grep(pattern = 'secondary', x = file_lines, value = TRUE, fixed = TRUE)[1]
    # mapped_secondary <- unlist(strsplit(x = mapped_secondary_line, split = ' '))[1]
    # mapped_total <- as.numeric(mapped_all) - as.numeric(mapped_secondary)
    
    # mapped to two different chromosomes
    chimeric_line <- grep(pattern = 'mate mapped to a different chr', x = file_lines, value = TRUE, fixed = TRUE)[1]
    chimeric <- as.numeric(unlist(strsplit(x = chimeric_line, split = ' '))[1])
    
    # the read is mapped in a proper pair
    properly_paired_line <- grep(pattern = " properly paired (", x = file_lines, value = TRUE, fixed = TRUE)[1]
    properly_paired <- as.numeric(unlist(strsplit(x = properly_paired_line, split = ' '))[1])
    
    
    return(list(
        SequencedPairedReads = total_paired,
        TotalMappedReads = mapped_all,
        ProperlyPairedReads = properly_paired,
        ProperlyPairedPcnt = properly_paired / total_paired,
        ChimericReads = chimeric
    ))
}

# ~~~~~~ RUN ~~~~~~ #
args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]

flagstat_df <- data.frame(read.flagstat(input_file))

write.table(x = flagstat_df, file = output_file, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# save.image(file = 'loaded.Rdata')
