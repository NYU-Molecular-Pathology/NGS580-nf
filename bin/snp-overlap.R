#!/usr/bin/env Rscript
# Script for overlapping homozygous SNPs from two annotation tables
library("ggplot2")
args <- commandArgs(TRUE)
tumor_tsv <- args[1]
normal_tsv <- args[2]
output_matrix <- args[3]
output_table <- args[4]
output_aggr_table <- args[5]
output_plot <- args[6]

tumor_df <- read.delim(tumor_tsv, header = T, sep = "\t", na.strings = '.', stringsAsFactors = FALSE)
normal_df <- read.delim(normal_tsv, header = T, sep = "\t", na.strings = '.', stringsAsFactors = FALSE)

save.image("loaded.Rdata")

# add a unique identifier to each variant
tumor_df[['VariantID']] <- sprintf("%s:%s:%s>%s",tumor_df[['CHROM']],tumor_df[['POS']],tumor_df[['REF']],tumor_df[['ALT']])
normal_df[['VariantID']] <- sprintf("%s:%s:%s>%s",normal_df[['CHROM']],normal_df[['POS']],normal_df[['REF']],normal_df[['ALT']])

save.image("loaded.Rdata")

# make sure that all VariantIDs are unique
if( any(duplicated(tumor_df[['VariantID']])) ) {
    print("ERROR: Duplicated VariantID's;")
    print(which(duplicated(tumor_df[['VariantID']])))
    quit(status = 1)
}
if( any(duplicated(normal_df[['VariantID']])) ) {
    print("ERROR: Duplicated VariantID's;")
    print(which(duplicated(normal_df[['VariantID']])))
    quit(status = 1)
}

## Functions to get elements from Venndiagram object in R
Intersect <- function (x) {
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's.
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

# make a list of the variants for overlapping
variant_list <- list(tumor=tumor_df[['VariantID']],normal=normal_df[['VariantID']])

# Gets all combinations of elements in a list
combs <- unlist(lapply(X = 1:length(variant_list), 
                       FUN = function(j) combn(names(variant_list), j, simplify = FALSE)), recursive = FALSE)

# get the elements that belong to each combination of groupings
elements <- lapply(X = combs, 
                   FUN = function(i) Setdiff(variant_list[i], variant_list[setdiff(names(variant_list), i)]))

# get the number of elements in each grouping
n.obs <- sapply(elements, length)

# generate a sequence of numbers up to the length of the largest set of elements
seq.max <- seq_len(max(n.obs))

# generate a matrix with the elements in separate columns
mat <- sapply(elements, "[", i = seq.max)

# generate some cleaner column names for the matrix
new_names <- sapply(combs, function(x){
    if ( length(x) > 1){
        name <- paste(x, collapse = '.')
        label <- sprintf("overlap: %s", name)
        return( label )
    } else{
        return(x)
    }
})
colnames(mat) <- new_names

write.table(x = mat, file = output_matrix, quote = F,sep = "\t", row.names = F)

# convert the matrix to a long format dataframe
# NOTE: this can get slow for large numbers of combinations...
overlap_df <- data.frame()
for( i in seq(length(colnames(mat))) ){
    # make dataframe
    df <- as.data.frame(elements[[i]])

        # make combination label column
    df[["comb"]] <- colnames(mat)[i]
    
    # append to full dataframe
    if(nrow(overlap_df) < 1){
      overlap_df <- df
    } else {
      overlap_df <- rbind(overlap_df, df)
    }
}

# rename the first column holding the Variant IDs
names(overlap_df)[1] <- "VariantID"

# add dummy variable for aggregating
overlap_df[['n']] <- 1

# add dummy variable for plotting
overlap_df[['all']] <- '.'

write.table(x = overlap_df, file = output_table, quote = FALSE, sep = '\t', row.names = FALSE)

# get total number of variants
total_num_variants <- length(unique(as.character(overlap_df[["VariantID"]])))

# aggregate on the percent of entries in each grouping
overlap_aggr <- aggregate(n ~ comb, data = overlap_df, FUN = sum)
overlap_pcnt <- aggregate(n ~ comb, data = overlap_df, FUN = function(x){
    return( round( (sum(x) / total_num_variants) * 100, digits = 1) )
})
names(overlap_pcnt) <- c("comb", "pcnt")
overlap_aggr <- merge(overlap_aggr, overlap_pcnt)
overlap_aggr[["all"]] <- '.'

write.table(x = overlap_aggr, file = output_aggr_table, quote = FALSE, sep = '\t', row.names = FALSE)

pdf(file = output_plot)
ggplot(data = overlap_aggr, aes(x = all, y = pcnt, fill = comb)) +
    geom_bar(stat = "identity", position="stack") +
    theme_bw() +
    ggtitle("Homozygous SNP Overlap") +
    theme(
        panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()
          )
dev.off()

save.image("final.Rdata")
