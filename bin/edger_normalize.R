#!/usr/bin/env Rscript

library(edgeR)
library(optparse)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

get_args <- function() {

    option_list <- list(
        make_option("--count-file", dest = 'count_file', help = "Input file path")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Normalize counts using edgeR"
        ))

    return(args)
}

get_normalized_counts <- function(count_file) {
    count_data <- read.csv(args$count_file, row.names = 1)
    count_data_matrix <- as.matrix(count_data)

    # Add a small pseudocount of 0.01 to avoid zero counts
    count_data_matrix[count_data_matrix == 0] <- 0.01

    dge <- DGEList(counts = count_data_matrix)
    rownames(dge) <- rownames(count_data_matrix)
    colnames(dge) <- colnames(count_data_matrix)

    #normalization
    dge <- calcNormFactors(dge, method="TMM")
    normalized_counts <- cpm(dge)

    return(normalized_counts)
}

export_data <- function(normalized_counts) {
    normalized_filename <- sub(".csv", "_normalized.csv", basename(args$count_file))
    write.table(normalized_counts, normalized_filename,, row.names = TRUE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

normalized_counts <- get_normalized_counts(args$count_file)

export_data(normalized_counts)
