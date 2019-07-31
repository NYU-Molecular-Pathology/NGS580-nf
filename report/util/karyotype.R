#!/usr/bin/env Rscript
library("ggplot2")
library("scales")

get_chrom_sizes <- function(){
# chrom_sizes <- read.delim(file = "chrom_sizes.hg19.txt", header = TRUE, sep = '\t')
chrom_sizes <- structure(list(chromosome = structure(c(23L, 1L, 12L, 16L, 17L,
18L, 19L, 20L, 21L, 22L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L,
11L, 13L, 14L, 15L, 24L, 25L),
.Label = c("chr1", "chr10", "chr11",
"chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
"chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5",
"chr6", "chr7", "chr8", "chr9", "chrM", "chrX", "chrY"), class = "factor"),
size = c(16571L, 249250621L, 243199373L, 198022430L, 191154276L,
180915260L, 171115067L, 159138663L, 146364022L, 141213431L,
135534747L, 135006516L, 133851895L, 115169878L, 107349540L,
102531392L, 90354753L, 81195210L, 78077248L, 59128983L, 63025520L,
48129895L, 51304566L, 155270560L, 59373566L)),
class = "data.frame", row.names = c(NA,
-25L))
}

get_centromeres <- function(){
# centromeres <- read.delim(file = "centromeres.hg19.txt", header = TRUE, sep = '\t')
centromeres <- structure(list(bin = c(23L, 20L, 2L, 1L, 14L, 16L, 1L, 14L, 1L,
1L, 10L, 1L, 15L, 13L, 1L, 1L, 11L, 13L, 1L, 1L, 1L, 12L, 10L,
10L), chromosome = structure(c(1L, 12L, 16L, 17L, 18L, 19L, 20L,
21L, 22L, 23L, 24L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L,
13L, 14L, 15L),
.Label = c("chr1", "chr10", "chr11", "chr12",
"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
"chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6",
"chr7", "chr8", "chr9", "chrX", "chrY"), class = "factor"),
start = c(121535434L,
92326171L, 90504854L, 49660117L, 46405641L, 58830166L, 58054331L,
43838887L, 47367679L, 58632012L, 10104553L, 39254935L, 51644205L,
34856694L, 16000000L, 16000000L, 17000000L, 35335801L, 22263006L,
15460898L, 24681782L, 26369569L, 11288129L, 13000000L), end = c(124535434L,
95326171L, 93504854L, 52660117L, 49405641L, 61830166L, 61054331L,
46838887L, 50367679L, 61632012L, 13104553L, 42254935L, 54644205L,
37856694L, 19000000L, 19000000L, 20000000L, 38335801L, 25263006L,
18460898L, 27681782L, 29369569L, 14288129L, 16000000L),
ix = c(1270L,
770L, 784L, 447L, 452L, 628L, 564L, 376L, 411L, 583L, 105L, 341L,
447L, 304L, 3L, 3L, 3L, 354L, 192L, 125L, 410L, 275L, 22L, 3L
),
n = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "N", class = "factor"),
size = c(3000000L, 3000000L, 3000000L, 3000000L, 3000000L,
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L,
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L,
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L,
3000000L), type = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L), .Label = "centromere", class = "factor"),
bridge = structure(c(1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "no", class = "factor")),
class = "data.frame", row.names = c(NA,
-24L))
}

get_chrom_order <- function(){
# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                 "chr22", "chrX", "chrY", "chrM")
}

get_karyotype_data <- function(){
    # return a list of the data needed to make the karyotype
    chrom_sizes <- get_chrom_sizes()
    centromeres <- get_centromeres()
    chrom_order <- get_chrom_order()
    
    chrom_key <- setNames(object = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                                                  12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                                  21, 22, 23, 24, 25)), 
                          nm = chrom_order)
    chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
    
    # convert the chromosome column in each dataset to the ordered factor
    chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                          levels = chrom_order)
    centromeres[["chromosome"]] <- factor(x = centromeres[["chromosome"]], 
                                          levels = chrom_order)
    
    data <- list(
        chrom_sizes = chrom_sizes,
        centromeres = centromeres, 
        chrom_key = chrom_key, 
        chrom_order = chrom_order
    )
    
    return(data)
}

karyotype_plot <- function(){
    # generates the base karyotype plot
    karyotype_data <- get_karyotype_data()
    chrom_sizes <- karyotype_data[["chrom_sizes"]]
    chrom_key <- karyotype_data[["chrom_key"]]
    centromeres <- karyotype_data[["centromeres"]]

    p <- ggplot(data = chrom_sizes) + 
        # base rectangles for the chroms, with numeric value for each chrom on the x-axis
        geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                      xmax = as.numeric(chromosome) + 0.2, 
                      ymax = size, ymin = 0), 
                  colour="black", fill = "white") + 
        # rotate the plot 90 degrees
        coord_flip() +
        # black & white color theme 
        theme(axis.text.x = element_text(colour = "black"), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank()) + 
        # give the appearance of a discrete axis with chrom labels
        scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
        # add bands for centromeres
        geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.2, 
                                          xmax = as.numeric(chromosome) + 0.2, 
                                          ymax = end, ymin = start)) +
        scale_y_continuous(labels = comma) +
        ylab("region (bp)") 

    return(p)
    
    # add variants code example: 
    # geom_point(data = annotations, 
    #            aes(x = as.numeric(chromosome),
    #                y = Start, color = "red"))
    # geom_rect(data = sample_cns, aes(xmin = as.numeric(chromosome) - 0.1, 
    #                                  xmax = as.numeric(chromosome) + 0.1, 
    #                                  ymax = end, ymin = start, fill = CNA)) + 
    # scale_fill_manual(values = group.colors) +
    # ggtitle("Variant Karotype") +
    # supress scientific notation on the y-axis
    # theme(legend.position="none")
}
