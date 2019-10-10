#!/usr/bin/env Rscript

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

dataframe_difference <- function(larger_df, smaller_df){
    # https://stackoverflow.com/questions/17427916/r-selecting-all-rows-from-a-data-frame-that-dont-appear-in-another
    # dataframe_difference(merged_df, avinput)
    cols_to_use <- names(smaller_df)[which(names(larger_df) %in% names(smaller_df))]
    x <- rbind(larger_df[, cols_to_use], smaller_df)
    return(x[! duplicated(x, fromLast=TRUE) & seq(nrow(x)) <= nrow(larger_df), ])
}

# ~~~~~~ RUN ~~~~~~ #
args <- commandArgs(TRUE)
# vcf_tsv_file <- "NC-HAPMAP.updated.tsv"
# annovar_txt_file <- "NC-HAPMAP.hg19_multianno.txt"
# avinput_file <- 'NC-HAPMAP.avinput.reformat.tsv'
# output_file <- "output.tsv"
vcf_tsv_file <- args[1]
annovar_txt_file <- args[2]
avinput_file <- args[3]
output_file <- args[4]

save.image('args.Rdata')

message(">>> Loading variant tables to be merged")

annovar <- read.ANNOVAR.vcf_txt(file = annovar_txt_file)

save.image('loaded.Rdata')

# need to convert all nucleotides to upper case; 
# some masked nucleotides from .vcf can be lower case, 
# but get converted to upper case in vcf .tsv table by GATK
any_masked_nucleotides_annovar <- any(grepl(pattern = 't|c|a|g',
                                            x = c(annovar[["Ref"]],
                                                  annovar[["Alt"]]) ))

if(any_masked_nucleotides_annovar) {
    message(">>> Detected lower-case masked nucleotides in annovar, converting to upper-case nucleotides")
    annovar[["Ref"]] <- toupper(annovar[["Ref"]])
    annovar[["Alt"]] <- toupper(annovar[["Alt"]])
}

save.image('loaded.Rdata')

vcf_table <- read.delim(file = vcf_tsv_file, 
                        header = TRUE, 
                        sep = '\t', 
                        stringsAsFactors = FALSE, 
                        colClasses = "character")
save.image('loaded.Rdata')

avinput <- read.delim(file = avinput_file, 
                      header = TRUE, 
                      sep = '\t', 
                      stringsAsFactors = FALSE, 
                      colClasses = "character")
save.image('loaded.Rdata')

# need to convert all nucleotides to upper case; 
# some masked nucleotides from .vcf can be lower case, 
# but get converted to upper case in vcf .tsv table by GATK
any_masked_nucleotides_avinput <- any(grepl(pattern = 't|c|a|g',
                                    x = c(avinput[["REF"]],
                                          avinput[["ALT"]],
                                          avinput[["Ref"]],
                                          avinput[["Alt"]]) ))

if(any_masked_nucleotides_avinput) {
    message(">>> Detected lower-case masked nucleotides in avinput, converting to upper-case nucleotides")
    avinput[["REF"]] <- toupper(avinput[["REF"]])
    avinput[["ALT"]] <- toupper(avinput[["ALT"]])
    avinput[["Ref"]] <- toupper(avinput[["Ref"]])
    avinput[["Alt"]] <- toupper(avinput[["Alt"]])
}

save.image('loaded.Rdata')


message(sprintf(">>> annovar: Number of rows: %s", nrow(annovar)))
message(sprintf(">>> avinput: Number of rows: %s", nrow(avinput)))
message(sprintf(">>> vcf_table: Number of rows: %s", nrow(vcf_table)))

# check that all df's have the same number of rows
if(! all_equal(nrow(annovar), nrow(avinput), nrow(vcf_table)) ) {
    message(">>> ERROR: Loaded dataframes have unequal numbers of rows")
    quit(status = 1)
}

message(">>> Merging tables")

# merge tables
merged_df <- Reduce(function(x, y){ merge(x, y, all = TRUE) }, list(avinput, vcf_table, annovar))

save.image('merged.Rdata')

message(sprintf(">>> merged_df: Number of rows: %s", nrow(merged_df)))

# check that all df's have the same number of rows
if(! all_equal(nrow(annovar), nrow(avinput), nrow(vcf_table), nrow(merged_df)) ) {
    message(">>> ERROR: Dataframes have unequal numbers of rows after merging")
    message(">>> Discrepant rows:")
    print(dataframe_difference(merged_df, avinput))
    message(">>> Rows with possible errors:")
    print(merged_df[merged_df[["POS"]] %in% dataframe_difference(merged_df, avinput)[["POS"]], names(dataframe_difference(merged_df, avinput)) ])
    quit(status = 1)
}

write.table(x = merged_df, file = output_file, sep = '\t', quote = FALSE, na = '.', row.names = FALSE, col.names = TRUE)
