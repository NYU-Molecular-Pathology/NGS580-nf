#!/usr/bin/env Rscript

# script to merge the ANNOVAR annotation tables plus GATK VariantsToTable tsv file for Pindel output
# because Pindel output needs special handling due to discrepancies compared to other outputs

# ~~~~~~ CUSTOM FUNCTIONS ~~~~~~~ #
read.ANNOVAR.vcf_txt <- function(file, na_string = '.', drop_Otherinfo = TRUE){
    # read the '_multianno.txt' file produced by ANNOVAR table_annovar.pl using the '--vcfinput' arg
    # tab-delimited text file has more columns than column names
    # output dataframe with colnames present
    
    # read in the file without headers; headers in first row
    tmp_annovar_df <- read.delim(file = file, 
                                 header = FALSE, 
                                 sep = '\t', 
                                 check.names = FALSE, 
                                 fill = TRUE, stringsAsFactors = FALSE, 
                                 na.strings = na_string)
    # get the colnames that are present
    av_colnames <- tmp_annovar_df[1,][which(! is.na(tmp_annovar_df[1,]) & ! tmp_annovar_df[1,] == "" )]
    # create new df without colnames
    tmp_annovar_df2 <- tmp_annovar_df[2:nrow(tmp_annovar_df),]
    # add the colnames
    names(tmp_annovar_df2)[1:length(av_colnames)] <- av_colnames
    
    # drop columns after 'Otherinfo', if its present
    if(isTRUE(drop_Otherinfo)){
        if (any(grepl(pattern = 'Otherinfo', x = names(tmp_annovar_df2)))){
            Otherinfo_index <- which(names(tmp_annovar_df2) == 'Otherinfo')
            tmp_annovar_df2 <- tmp_annovar_df2[, names(tmp_annovar_df2)[1:Otherinfo_index - 1]]
        }
    }
    
    return(tmp_annovar_df2)
}

all_equal <- function(...){
    # returns TRUE or FALSE if all the elements passed are equal
    arguments <- list(...)
    TF <- vapply(1:(length(arguments)-1),
                 function(n) identical(arguments[[n]], arguments[[n+1]]),
                 logical(1))
    return(all(TF))
}

# ~~~~~~ RUN ~~~~~~ #
args <- commandArgs(TRUE)
# vcf_tsv_file <- "NC-HAPMAP.updated.tsv"
# annovar_txt_file <- "NC-HAPMAP.hg19_multianno.txt"
# avinput_file <- 'NC-HAPMAP.avinput.reformat.tsv'
# output_file <- "output.tsv"
vcf_file <- args[1] # the original sample .vcf
vcf_tsv_file <- args[2] # tsv version of the .vcf
annovar_txt_file <- args[3] # ANNOVAR annotation output
avinput_file <- args[4] # ANNOVAR avinput file with only selected columns for merge
output_file <- args[5] # name of output final merged table file

save.image('args.Rdata')

message(">>> Loading variant tables to be merged")

annovar <- read.ANNOVAR.vcf_txt(file = annovar_txt_file)
save.image('loaded.Rdata')

vcf_table <- read.delim(file = vcf_tsv_file, 
                        na.strings = '.',
                        header = TRUE, 
                        sep = '\t', 
                        stringsAsFactors = FALSE, 
                        colClasses = "character")
save.image('loaded.Rdata')

# read in the raw .vcf file
# skip header rows; hard code them in later cause we already know what they should be
vcf_raw_df <- read.delim(file = vcf_file, 
                         header = FALSE, 
                         na.strings = '.',
                         sep = '\t', 
                         comment.char = '#', 
                         colClasses = "character",
                         col.names = c("CHROM", "POS", "ID", "REF", "ALT", 
                                       "QUAL", "FILTER", "INFO","FORMAT", "NORMAL", "TUMOR"))

avinput <- read.delim(file = avinput_file, 
                      header = FALSE, 
                      sep = '\t', 
                      stringsAsFactors = FALSE, 
                      colClasses = "character",
                      na.strings = '.',
                      col.names = c("Chr", "Start", "End", "Ref", "Alt", # from ANNOVAR
                                    "CHROM", "POS", "ID", "REF", "ALT", # from vcf
                                    "QUAL", "FILTER", "INFO","FORMAT", "NORMAL", "TUMOR"))
save.image('loaded.Rdata')

message(sprintf(">>> annovar: Number of rows: %s", nrow(annovar)))
message(sprintf(">>> avinput: Number of rows: %s", nrow(avinput)))
message(sprintf(">>> vcf_table: Number of rows: %s", nrow(vcf_table)))
message(sprintf(">>> vcf_raw_df: Number of rows: %s", nrow(vcf_raw_df)))


# check that all df's have the same number of rows
if(! all_equal(nrow(annovar), nrow(avinput), nrow(vcf_table), nrow(vcf_raw_df)) ) {
    print("ERROR: Loaded dataframes have unequal numbers of rows")
    quit(status = 1)
}

message(">>> Merging tables")

# # METHODS FOR IDENTIFYING DUPLICATES:
# # NEED THESE NOTES FOR DEBUGGING
# y <- rbind(vcf_raw_df[, intersect(names(vcf_table), names(vcf_raw_df))], vcf_table[, intersect(names(vcf_table), names(vcf_raw_df)) ])
# y[! duplicated(y, fromLast=TRUE) & seq(nrow(y)) <= nrow(vcf_table), ]

# x <- merge(vcf_table, vcf_raw_df, all = TRUE)
# nrow(x)
# x[duplicated(x[, names(vcf_raw_df)]), ]
# x <- merge(x, avinput, all = TRUE)
# nrow(x)
# x <- merge(x, annovar, all = TRUE)
# nrow(x)
# x[duplicated(x[, names(annovar)]), ]

# x <- merge(x, vcf_table, all = TRUE)
# nrow(x)
# z <- rbind(x[, intersect(names(vcf_raw_df), names(x))], vcf_raw_df[, intersect(names(vcf_raw_df), names(x)) ])
# nrow(z)
# z[! duplicated(z, fromLast=TRUE) & seq(nrow(z)) <= nrow(vcf_raw_df), ]

# cols <- intersect(names(avinput), names(vcf_table))
# x <- merge(vcf_table, avinput, all = FALSE)
# x <- rbind(avinput[, which(names(avinput) %in% names(vcf_table))], vcf_table[, which(names(vcf_table) %in% names(avinput))])
# x[! duplicated(x, fromLast=TRUE) & seq(nrow(x)) <= nrow(vcf_table), ]

# merge tables
merged_df <- Reduce(function(x, y){ merge(x, y, all = TRUE) }, list(vcf_raw_df, vcf_table, avinput, annovar))

save.image('merged.Rdata')

message(sprintf(">>> merged_df: Number of rows: %s", nrow(merged_df)))

# check that all df's have the same number of rows
if(! all_equal(nrow(annovar), nrow(avinput), nrow(vcf_table), nrow(merged_df)) ) {
    print("ERROR: Dataframes have unequal numbers of rows after merging")
    quit(status = 1)
}

keep_cols <- names(merged_df)[which(!names(merged_df) %in% c('INFO', 'FORMAT', 'NORMAL', 'TUMOR'))]
write.table(x = merged_df[, keep_cols], file = output_file, sep = '\t', quote = FALSE, na = '.', row.names = FALSE, col.names = TRUE)

