#!/usr/bin/env Rscript

library(DESeq2)
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
        description = "Normalize counts using DESeq2"
        ))

    return(args)
}

get_normalized_counts <- function(count_file) {
    print(count_file)
    count_data <- read.csv(args$count_file, row.names = 1)

    count_data_matrix <- as.matrix(count_data)
    colData <- DataFrame(row.names = colnames(count_data_matrix))

    # Add a small pseudocount of 1 (it has to be integer...) to avoid zero counts
    count_data_matrix[count_data_matrix == 0] <- 1

    dds <- DESeqDataSetFromMatrix(countData = count_data_matrix, colData = colData, design = ~ 1)

    # perform normalization
    dds <- estimateSizeFactors(dds)

    normalized_counts <- counts(dds, normalized = TRUE)

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
