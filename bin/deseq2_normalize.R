#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

library(DESeq2)
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
        description = "Normalize counts using DESeq2"
        ))

    return(args)
}

prefilter_counts <- function(dds, design_data) {
    # see https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
    # getting size of smallest group
    group_sizes <- table(design_data$condition)
    smallest_group_size <- min(group_sizes)
    # keep genes with at least 10 counts over a certain number of samples
    keep <- rowSums(counts(dds) >= 10) >= smallest_group_size
    dds <- dds[keep,]
    return(dds)
}

get_normalized_counts <- function(dds) {
    # perform normalization
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
    return(normalized_counts)
}

filter_out_genes_with_at_least_one_zero <- function(count_matrix, normalized_counts) {
    # get gene IDs corresponding to rows with only non-zero counts
    non_zero_rows <- rownames(count_matrix[apply(count_matrix!=0, 1, all),])
    # filter out genes with zero counts
    normalized_counts <- normalized_counts[rownames(normalized_counts) %in% non_zero_rows, ]
    return(normalized_counts)
}

get_log2_cpm_counts <- function(normalized_counts, count_data) {
    # calculate total counts per sample (library size)
    library_sizes <- colSums(count_data)
    # convert normalized counts to CPM
    cpm_counts <- t(t(normalized_counts) / library_sizes * 1e6)
    cpm_counts <- log2(cpm_counts + 1)
    return(cpm_counts)
}

get_normalized_cpm_counts <- function(count_file, design_file) {

    print(paste('Normalizing counts in:', count_file))

    count_data <- read.csv(count_file, row.names = 1)
    design_data <- read.csv(design_file)

    design_data <- design_data[design_data$sample %in% colnames(count_data), ]

    count_matrix <- as.matrix(count_data)

    col_data <- data.frame(
        row.names = design_data$sample,
        condition = factor(design_data$condition)
    )

    # check if the column names of count_matrix match the row names of col_data
    if (!all(colnames(count_matrix) == rownames(col_data))) {
        stop("Sample names in the count matrix do not match the design data.")
    }

    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)

    # pre-filter genes with low counts
    # not absolutely necessary, but good practice to avoid RAM issues
    dds <- prefilter_counts(dds, design_data)

    normalized_counts <- get_normalized_counts(dds)
    #print(length(rownames(normalized_counts)))

    normalized_counts <- filter_out_genes_with_at_least_one_zero(count_matrix, normalized_counts)
    #print(length(rownames(normalized_counts)))

    cpm_counts <- get_log2_cpm_counts(normalized_counts, count_data)

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
