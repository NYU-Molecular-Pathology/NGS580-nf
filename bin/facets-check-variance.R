#!/usr/bin/env Rscript
# Script to check the variance on a snp-pileup
# as pre-filter step before running full Facets analysis

library(facets)

args <- commandArgs(TRUE)
snp_pileup_file <- args[1] # standard snp-pilup output
output_var_file <- args[2]
save.image("args.Rdata")

datafile = snp_pileup_file # path to the snp-pileup file
save.image("loaded.Rdata")

rcmat = readSnpMatrix(datafile) #read in the matrix
save.image("loaded.Rdata")

preproc_sample = preProcSample(rcmat)
save.image("loaded.Rdata")

oo=procSample(preproc_sample,cval=150) ##  Lower cval lead to higher sensitivity for small changes
save.image("loaded.Rdata")

snp_var <- var(oo$jointseg$cnlr,na.rm=T)

fileConn<-file(output_var_file)
writeLines(as.character(snp_var), fileConn)
close(fileConn)

