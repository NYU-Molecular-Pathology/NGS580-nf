#!/usr/bin/env Rscript
library("ggplot2")
args <- commandArgs(TRUE)
tumor_tsv <- args[1]
normal_tsv <- args[2]

tumor_df <- read.delim(tumor_tsv, header = T, sep = "\t", na.strings = '.')
normal_df <- read.delim(normal_tsv, header = T, sep = "\t", na.strings = '.')

tumor_df[['VariantID']] <- sprintf("%s%s%s%s",tumor_df[['CHROM']],tumor_df[['POS']],tumor_df[['REF']],tumor_df[['ALT']])
normal_df[['VariantID']] <- sprintf("%s%s%s%s",normal_df[['CHROM']],normal_df[['POS']],normal_df[['REF']],normal_df[['ALT']])

mincov <- 200 #DP
minvaf <- 0.98 #FREQ

tumor_df <- tumor_df[which(tumor_df[['DP']] > mincov), ]
tumor_df <- tumor_df[which(tumor_df[['FREQ']] > minvaf), ]

normal_df <- normal_df[which(normal_df[['DP']] > mincov), ]
normal_df <- normal_df[which(normal_df[['FREQ']] > minvaf), ]

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


variant_list <- list(tumor=tumor_df[['VariantID']],normal=normal_df[['VariantID']])
## Gets all combinations of elements in a list
combs <- unlist(lapply(1:length(variant_list), function(j) combn(names(variant_list), j, simplify = FALSE)), recursive = FALSE)
elements <- lapply(combs, function(i) Setdiff(variant_list[i], variant_list[setdiff(names(variant_list), i)]))
n.obs <- sapply(elements, length)
seq.max <- seq_len(max(n.obs))
mat <- sapply(elements, "[", i = seq.max)
colnames(mat) <- combs
#colnames(mat) <- gsub(",",".vs.",colnames(mat))
# write.table(mat,file="combinations.tsv", quote=F,sep="\t",row.names=F)


overlap_df <- data.frame()
for( i in seq(length(combs)) ){
  print(class(elements[[i]]))
  df <- as.data.frame(elements[[i]])
  df[["comb"]] <- paste(combs[[i]], collapse = '.')
  
  if(nrow(overlap_df) < 1){
    overlap_df <- df
  } else {
    overlap_df <- rbind(overlap_df, df)
  }
}
names(overlap_df)[1] <- "VariantID"
# add dummy variable for plotting
overlap_df[['all']] <- '.'

ggplot(data = overlap_df, aes(x = all, fill = comb)) + geom_bar(position="stack")
save.image("data.Rdata")

