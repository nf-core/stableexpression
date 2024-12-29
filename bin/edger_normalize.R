#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

library(edgeR)
library(optparse)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

get_args <- function() {

    option_list <- list(
        make_option("--counts", dest = 'count_file', help = "Path to input count file"),
        make_option("--design", dest = 'design_file', help = "Path to input design file")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Normalize counts using edgeR"
        ))

    return(args)
}

replace_zero_counts_with_pseudocounts <- function(count_data_matrix) {
    # Add a small pseudocount of 0.01 to avoid zero counts
    count_data_matrix[count_data_matrix == 0] <- 0.01
    return(count_data_matrix)
}

filter_out_lowly_expressed_genes <- function(dge) {
    # filter the dataframe to exclude rows where the mean is 0
    # filter out lowly expressed genes
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    return(dge)
}

get_non_zero_rows <- function(count_matrix) {
    # get gene IDs corresponding to rows with only non-zero counts
    non_zero_rows <- rownames(count_matrix[apply(count_matrix!=0, 1, all),])
    return(non_zero_rows)
}

filter_out_genes_with_at_least_one_zero <- function(non_zero_rows, cpm_counts) {
    # filter out genes with zero counts
    cpm_counts <- cpm_counts[rownames(cpm_counts) %in% non_zero_rows, ]
    return(cpm_counts)
}

get_log2_cpm_counts <- function(dge) {
    cpm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
    return(cpm_counts)
}


get_normalized_cpm_counts <- function(count_file, design_file) {

    print(paste('Normalizing counts in:', count_file))

    count_data <- read.csv(args$count_file, row.names = 1)
    design_data <- read.csv(design_file)

    design_data <- design_data[design_data$sample %in% colnames(count_data), ]
    group <- factor(design_data$condition)
    count_matrix <- as.matrix(count_data)

    non_zero_rows <- get_non_zero_rows(count_matrix)

    count_matrix <- replace_zero_counts_with_pseudocounts(count_matrix)

    dge <- DGEList(counts = count_matrix, group = group)
    rownames(dge) <- rownames(count_matrix)
    colnames(dge) <- colnames(count_matrix)

    dge <- filter_out_lowly_expressed_genes(dge)

    # normalization
    dge <- calcNormFactors(dge, method="TMM")

    cpm_counts <- get_log2_cpm_counts(dge)

    cpm_counts <- filter_out_genes_with_at_least_one_zero(non_zero_rows, cpm_counts)

    return(cpm_counts)
}

export_data <- function(cpm_counts, filename) {
    filename <- sub("\\.csv$", ".log_cpm.csv", filename)
    print(paste('Exporting normalized counts per million to:', filename))
    write.table(cpm_counts, filename, sep = ',', row.names = TRUE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

cpm_counts <- get_normalized_cpm_counts(args$count_file, args$design_file)

export_data(cpm_counts, basename(args$count_file))
