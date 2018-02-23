#!/usr/bin/env Rscript
args <- commandArgs(T)

sampleID <- args[1]
vcf_file <- args[2]
# vcf_file <- "pipeline_output/VCF-GATK-HC/SC-SERACARE.vcf"


message(sprintf("vcf_file: %s", vcf_file))
message(sprintf("sampleID: %s", sampleID))

save.image()

library("BSgenome.Hsapiens.UCSC.hg19")
library("deconstructSigs")

vcf_colnames <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sampleID)
variants <- read.delim(file = vcf_file, header = FALSE, sep = '\t',
                       comment.char = '#', col.names = vcf_colnames, check.names = FALSE)

# add sample ID column
variants[["Sample"]] <- rep(sampleID, nrow(variants))

# keep only entries with chroms in the reference data
variants <- variants[which(as.character(variants[["CHROM"]]) %in% seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)), ]


# need at least 55 variants per sample
if (nrow(variants) < 55) {
    message(sprintf('There are fewer than 55 variants for sample %s', sampleID))
    quit(status = 11)
}

# convert to signatures format
sigs.input <- mut.to.sigs.input(mut.ref = variants,
                                sample.id = "Sample",
                                chr = "CHROM",
                                pos = "POS",
                                ref = "REF",
                                alt = "ALT")

# make the signatures
signatures <- whichSignatures(tumor.ref = sigs.input,
                              signatures.ref = signatures.cosmic,
                              sample.id = sampleID,
                              contexts.needed = TRUE,
                              tri.counts.method = 'default')

# save signatures
signatures_output_file <- sprintf(sprintf('%s_signatures.Rds', sampleID))
saveRDS(object = signatures, file = signatures_output_file, compress = TRUE)

# make plots
signatures_plot_file <- sprintf('%s_signatures.pdf', sampleID)
pdf(file = signatures_plot_file)
print(plotSignatures(signatures, sub = 'signatures.cosmic'))
dev.off()

signatures_pie_plot_file <- sprintf('%s_signatures_pie.pdf', sampleID)
pdf(file = signatures_pie_plot_file)
print(makePie(signatures, sub = 'signatures.cosmic'))
dev.off()
