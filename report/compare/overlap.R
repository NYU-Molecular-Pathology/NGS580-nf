#!/usr/bin/env Rscript
# Module for overlap functions, used to overlap lists of variants

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

create_overlaps <- function(variant_list, col_name = "VariantID"){
    # convert a list of character vectors into a dataframe of overlap combinations
    # variant_list <- list(tumor=tumor_df[['VariantID']],normal=normal_df[['VariantID']])
    
    # make sure there are no duplicates
    for( name in names(variant_list)){
        if( any(duplicated(variant_list[[name]])) ) {
            stop("Duplicate labels present in variant list")
        }
    }
    
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

    # convert the matrix to a long format dataframe
    # NOTE: this can get slow for large numbers of combinations...
    overlap_df <- data.frame()
    for( i in seq(length(colnames(mat))) ){
        # make dataframe
        # # check that there are elements to make df from...
        if (length(elements[[i]]) < 1){
            # create dataframe with no rows
            df <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), 
                           c(col_name, colnames(mat)[i]))
        } else {
            df <- as.data.frame(elements[[i]])
            # make combination label column
            df[["comb"]] <- colnames(mat)[i]
            names(df)[1] <- col_name
        }
        
        # append to full dataframe
        if(nrow(overlap_df) < 1){
            overlap_df <- df
        } else {
            overlap_df <- rbind(overlap_df, df)
        }
    }
    
    return(overlap_df)
}

aggregate_overlaps <- function(overlap_df, col_name = "VariantID"){
    # create an aggregate table of the overlaps, for plotting
    # add dummy variable for aggregating
    overlap_df[["n"]] <- 1
    
    # add dummy variable for plotting
    overlap_df[["all"]] <- '.'
    
    # get total number of variants
    total_num_unpaired_variants <- length(unique(as.character(overlap_df[[col_name]])))
    
    # # aggregate on the percent of entries in each grouping
    overlap_aggr <- aggregate(n ~ comb, 
                              data = overlap_df,
                              FUN = sum)
    
    overlap_pcnt <- aggregate(n ~ comb, 
                              data = overlap_df, 
                              FUN = function(x){
                                  return( round( (sum(x) / total_num_unpaired_variants) * 100, digits = 1) )
                              })
    names(overlap_pcnt) <- c("comb", "pcnt")
    overlap_aggr <- merge(overlap_aggr, overlap_pcnt)
    overlap_aggr[["all"]] <- '.'
    return(overlap_aggr)
}

# 
# pdf(file = output_plot)
# ggplot(data = overlap_aggr, aes(x = all, y = pcnt, fill = comb)) +
#     geom_bar(stat = "identity", position="stack") +
#     theme_bw() +
#     ggtitle("Homozygous SNP Overlap") +
#     theme(
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank()
#     )
# dev.off()
# 
# save.image("final.Rdata")
