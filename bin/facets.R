#!/usr/bin/env Rscript

library(facets)

args <- commandArgs(TRUE)
snp_pileup_file <- args[1] # standard snp-pilup output
output_pdf_file <- args[2] # file to save Plotly PDF to
output_seg_file <- args[3] # file to save the csv CNV segment
save.image("args.Rdata")

datafile = snp_pileup_file # path to the snp-pileup file
save.image("loaded.Rdata")

rcmat = readSnpMatrix(datafile) #read in the matrix
save.image("loaded.Rdata")

preproc_sample = preProcSample(rcmat)
save.image("loaded.Rdata")
## A bivariate genome segmentation is performed on logR and logOR by extending the CBS algotithm (Olshen et al., 2004; Venkatraman and Olshen, 2007) to the bivariate scenario using a Tsquared statistic for identifiying change points ##

oo=procSample(preproc_sample,cval=150) ##  Lower cval lead to higher sensitivity for small changes
save.image("loaded.Rdata")
## Call allele-specific copy number and associated cellular fraction, estimate tumor purity and ploidy

fit=emcncf(oo)
save.image("loaded.Rdata")
#head(fit$cncf) ## The segmentation result and the EM fit output looks like this
#fit$purity ## estimated sample purity
#fit$ploidy ## estimated sample ploidy

save.image("loaded.Rdata")
write.csv(fit$cncf,file=output_seg_file,quote = F,row.names = F)

pdf(output_pdf_file,width = 10, height = 10)
sname<- sprintf('ploidy= %.2f; purity= %.2f', fit$ploidy, fit$purity)
plotSample(x=oo,emfit=fit,sname=sname)
dev.off()

## The top panel of the figure displays logR with chromosomes alternating in blue and gray. The green line indicates the median logR in the sample. The purple line indicates the logR of the diploid state. The second panel displays logOR. Segment means are ploted in red lines. The third panel plots the total (black) and minor (red) copy number for each segment. The bottom bar shows the associated cellular fraction (cf). Dark blue indicates high cf. Light blue indicates low cf. Beige indicates a normal segment (total=2,minor=1) ##