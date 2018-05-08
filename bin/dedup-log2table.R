#!/usr/bin/env Rscript

# this script will convert the sambamaba dedup stderr messages into a table
# add more values to output columns as needed

# ~~~~~~ FUNCTIONS ~~~~~~~ #
read.dedup.log <- function(file){
    # reads the stderr message from sambamba dedup
    # file <- input_file
    file_lines <- readLines(file)
    dup_line <- grep(pattern = 'found.*duplicates', x = file_lines, value = TRUE)[1]
    dups <- as.numeric(gsub(pattern = '.*found ([0-9]*) duplicates.*', replacement = '\\1', x = dup_line))
    
    total_pairs_line <- grep(pattern = 'sorted.*end pairs', x = file_lines, value = TRUE)[1]
    total_pairs <- as.numeric(gsub(pattern = '.* ([0-9]*) .*', replacement = '\\1', x = total_pairs_line, perl = TRUE))
    
    return(list(
        InputReadPairs = total_pairs,
        Duplicates = dups
    ))
}


# ~~~~~~ RUN ~~~~~~ #
args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]

dedup_df <- data.frame(read.dedup.log(input_file))

write.table(x = dedup_df, file = output_file, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# save.image("loaded.Rdata")
