#!/usr/bin/env Rscript
# Create a Plotly plot for CNVkit output
# Authors: Stephen Kelly, Varshini Vasudevaraja
library("ggplot2")
library("plotly")

args <- commandArgs(TRUE)
cns_file <- args[1] # standard CNVkit output
output_html_file <- args[2] # file to save Plotly HTML to
output_pdf_file <- args[3] # file to save Plotly PDF to

cns <- read.table(cns_file,header = T,sep = "\t")

save.image("loaded.Rdata")

# # create ordered factor level for chroms
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                 "chr22", "chrX", "chrY", "chrM")
chrom_key <- setNames(object = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                              12, 13, 14, 15, 16, 17, 18, 19, 20,
                                              21, 22, 23, 24, 25)),
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
cns[["chromosome"]] <- factor(x = cns[["chromosome"]],
                                levels = chrom_order)
cns[["gain_loss"]] <- factor(ifelse(cns[["log2"]] > 0, "gain", "loss"))



cnvplot <- ggplot(data = cns,
                  aes(x = start,
                      xend = end,
                      y = log2,
                      yend = log2,
                      label = gene,
                      color = gain_loss
                  )) +
    geom_point() +
    # geom_segment(size = 2, lineend = "round") + # does not work right in ggplotly
    geom_hline(yintercept = 0) +
    facet_grid(.~chromosome, scales = "free_x", drop = FALSE) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid = element_blank()) +
    scale_y_continuous(limits = c( min(cns[["log2"]]) - 0.5, max(cns[["log2"]]) + 0.5)) +
    theme(legend.position="none") +
    scale_color_manual(values = c("green", "red"))


cnvplotly <- ggplotly(cnvplot, tooltip = "gene")

htmlwidgets::saveWidget(as_widget(cnvplotly), file = output_html_file, selfcontained = TRUE)

pdf(file = output_pdf_file, width = 15, height = 8)
print(cnvplot)
dev.off()

save.image("final.Rdata")