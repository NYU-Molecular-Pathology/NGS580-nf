#!/usr/bin/env Rscript
# USAGE: ./VAF-distribution-plot.R foo /ifs/data/molecpathlab/PNET_GYN/sns_WES/VCF-GATK-HC-annot/256.combined.txt
args <- commandArgs(T)

sample_ID <- args[1]
vcf_annot_file <- args[2]

library("ggplot2")
# ~~~~~ FUNCTIONS ~~~~~ #
filter_annotations <- function(df, dbSNP_colname = "dbSNP_147", NA_str = '.'){
    # chroms to subset by
    # keep_chroms <- c("chrY", "chrX", "chr1", "chr2")
    # subset for chroms
    # df <- df[which(df[["CHR"]] %in% keep_chroms), ]
    
    # keep only dbSNP entries
    df <- df[ which(df[[dbSNP_colname]] != NA_str), ]
    
    return(df)
}


# ~~~~~ LOAD DATA ~~~~ #
annot_df <- read.delim(file = vcf_annot_file, header = TRUE, sep = '\t', check.names = FALSE)

annot_df <- filter_annotations(annot_df)

# quit if no variants
if (nrow(annot_df) < 1) {
    message(sprintf('There are not enough variants for sample %s', sample_ID))
    quit(status = 11)
}

pdf_output_file <- sprintf("%s_vaf_dist.pdf", sample_ID)
plot_title <- sprintf('%s - Variant Frequency Distribution', sample_ID)

vaf_plot <- ggplot(data = annot_df, aes(x = FREQ, fill = CHR)) +
    geom_vline(xintercept = 0.333, color = "red", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0.666, color = "red", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0.5, color = "blue", linetype = "dashed", alpha = 0.5) +
    geom_histogram(alpha = 0.8, bins = 100) + 
    ggtitle(plot_title) +
    coord_cartesian(ylim = c(0, 250), xlim = c(0, 1)) + #  
    theme_bw() +
    facet_grid(CHR ~ .) +
    scale_x_continuous(breaks = seq(0, 1, .1)) 

pdf(file = pdf_output_file, width = 10, height = 60)
print(vaf_plot)
dev.off()

save.image()
