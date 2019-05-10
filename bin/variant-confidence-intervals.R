#!/usr/bin/env Rscript
# Script to calculate confidence intervals for the SeraCare variants and add the values to the table
library("binom")
args <- commandArgs(TRUE)
input_tsv <- args[1]
output_tsv <- args[2]
error_rate <- as.numeric(args[3]) # 0.02
power <- 0.95
alpha <- 0.05
conf_level <- 1 - alpha

variants <- read.delim(input_tsv, 
                       header = T, 
                       sep = "\t", 
                       na.strings = '.', 
                       stringsAsFactors = FALSE, 
                       check.names = FALSE)
save.image("loaded.Rdata")

# calculate the required coverage for accurate variant calling
variants[["SeraCare.ReqCov"]] <- apply(X = variants, MARGIN = 1, FUN = function(row){
    True_AF <- as.numeric(as.character(row["SeraCare.AF"]))
    
    coverage_required <- cloglog.sample.size(p.alt = True_AF, 
                                             p = error_rate, 
                                             power = power, 
                                             alpha = alpha)[["n"]]
    return(coverage_required)
})


# add the upper and lower Confidence Intervals
variants[["SeraCare.CI.Lower"]] <- apply(X = variants, MARGIN = 1, FUN = function(row){
    coverage <- as.numeric(as.character(row["DP"]))
    True_AF <- as.numeric(as.character(row["SeraCare.AF"]))
    
    intervals <- binom.confint(x = True_AF * coverage, 
                               n = coverage, 
                               conf.level = conf_level, 
                               methods = "cloglog")
    CI_lower <- intervals[["lower"]]
    return(CI_lower)
})

variants[["SeraCare.CI.Upper"]] <- apply(X = variants, MARGIN = 1, FUN = function(row){
    coverage <- as.numeric(as.character(row["DP"]))
    True_AF <- as.numeric(as.character(row["SeraCare.AF"]))
    
    intervals <- binom.confint(x = True_AF * coverage, 
                               n = coverage, 
                               conf.level = conf_level, 
                               methods = "cloglog")
    CI_upper <- intervals[["upper"]]
    return(CI_upper)
})

write.table(x = variants, 
            file = output_tsv, 
            sep = '\t', 
            quote = FALSE, 
            na = '.', 
            row.names = FALSE, 
            col.names = TRUE)

save.image("final.Rdata")
