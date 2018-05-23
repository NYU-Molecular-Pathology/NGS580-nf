#!/usr/bin/env Rscript

# used on Skynet skygpu10, requires "libxml2-devel" installed on the system
Sys.setenv(XML_CONFIG="/usr/bin/xml2-config") 

# base R packages
install.packages(c(
    "ggplot2"
    ), repos='http://cran.us.r-project.org', dependencies = TRUE)

# bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite(c(
    "BSgenome",
    "BSgenome.Hsapiens.UCSC.hg19",
    "GenomeInfoDb"
    ))

# base R packages that depend on bioconductor packages
install.packages(c(
"deconstructSigs"
), repos='http://cran.us.r-project.org', dependencies = TRUE)
