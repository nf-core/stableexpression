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
        make_option("--counts", dest = 'count_file', help = "Path to input count file")
        make_option("--design", dest = 'design_file', help = "Path to input design file")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Normalize counts using DESeq2"
        ))

    return(args)
}

get_normalized_counts <- function(count_file, design_file) {

    count_data <- read.csv(count_file, row.names = 1)
    design_data <- read.csv(design_file, row.names = 1)

    count_matrix <- as.matrix(count_data)

    # Ensure the row names of the colData match the sample names in the count matrix
    col_data <- data.frame(
        row.names = design_data$sample_names,  # Assuming 'sample_names' is the column name
        condition = factor(design_data$groups)  # Assuming 'groups' is the column name
    )

    # Check if the column names of count_matrix match the row names of col_data
    if (!all(colnames(count_matrix) == rownames(col_data))) {
        stop("Sample names in the count matrix do not match the design data.")
    }

    # Add a small pseudocount of 1 (it has to be integer...) to avoid zero counts
    count_matrix[count_matrix == 0] <- 1

    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)

    # Filter genes with low counts
    # Keep genes that have at least 10 counts in at least two samples
    dds_filtered <- dds[rowSums(counts(dds) >= 2) >= 2, ]

    # perform normalization
    dds <- estimateSizeFactors(dds)

    normalized_counts <- counts(dds, normalized = TRUE)

    return(normalized_counts)
}

export_data <- function(normalized_counts, file_stem) {
    normalized_filename <- sub(".csv", "_normalized.csv", file_stem)
    write.table(normalized_counts, normalized_filename,, row.names = TRUE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

normalized_counts <- get_normalized_counts(args$count_file)

export_data(normalized_counts, file_stem = basename(args$count_file))
