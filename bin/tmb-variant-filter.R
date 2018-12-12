#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
Rdata <- args[1]
output <- args[2]
anno_txt <- args[3] 

# values used in the ANNOVAR output for nastring
NA_vals <- c('.')  

# sampleID <- "RD-18-508-S12-24039.HaplotypeCaller"
# avinput <- "output/RD-18-508-S12-24039/RD-18-508-S12-24039.HaplotypeCaller.avinput"
# anno_txt <- "output/RD-18-508-S12-24039/RD-18-508-S12-24039.HaplotypeCaller.hg19_multianno.txt"

read.annot_txt <- function(file){
    # read ANNOVAR hg19_multianno.txt when 'Otherinfo' results in more columns than headers
    # get the first two lines
    con <- file(file,"r")
    first_lines <- readLines(con,n=2)
    close(con)
    
    # count number of columns in table
    num_cols <- length(unlist(strsplit(first_lines[[2]], '\t')))
    # get existing headers
    df_names <- unlist(strsplit(first_lines[[1]], '\t'))
    # number of existing names
    num_names <- length(df_names)
    # fill in missing names
    all_names <- c(df_names, sprintf('V%s', seq(from = num_names + 1, to = num_cols)))
    # load table starting with line 2, using new headers
    df <- read.delim(file = file, header = FALSE, sep = '\t', col.names = all_names, skip = 1, na.strings = NA_vals)
    return(df)
}

filter_variants <- function(df){
    # inspired by:
    # https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0424-2
    # https://doi.org/10.1186/s13073-017-0424-2
    # Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden. Chalmers ZR et. al.
    # Non-coding alterations were not counted. 
    # Alterations listed as known somatic alterations in COSMIC and truncations in tumor suppressor genes were not counted, since our assay genes are biased toward genes with functional mutations in cancer [63]. 
    # Alterations predicted to be germline by the somatic-germline-zygosity algorithm were not counted [64].
    # Alterations that were recurrently predicted to be germline in our cohort of clinical specimens were not counted. 
    # Known germline alterations in dbSNP were not counted. 
    # Germline alterations occurring with two or more counts in the ExAC database were not counted [65]. 
    
    # exclude non-coding entries
    df <- df[which(!is.na(df[["AAChange.refGene"]])), ]
    # df <- df[which(! df[["AAChange.refGene"]] %in% NA_vals ), ]
    # same as filtering for: 'exonic', 'splicing', 'exonic;splicing'# "Func.refGene"
    
    
    # exclude entries in COSMIC; keep only 'NA' entries
    df <- df[which(is.na(df[["cosmic70"]])), ] # [["cosmic70"]]

    # exclude entries in dbSNP; keep only 'NA' entries
    df <- df[which(is.na(df[["snp138"]])), ]
    
    # exclude entries in ExAC; keep only 'NA' or 0 entries
    df <- df[which(is.na(df[["ExAC_ALL"]]) | df[["ExAC_ALL"]] == 0 ), ]
    return(df)
}

anno_df_og <- read.annot_txt(anno_txt)
anno_df <- filter_variants(anno_df_og)

write.table(x = anno_df, file = output, sep = '\t', row.names = FALSE, col.names = TRUE, na = NA_vals[1])

save.image(Rdata)