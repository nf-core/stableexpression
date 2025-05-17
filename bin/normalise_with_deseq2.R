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

check_samples <- function(count_matrix, design_data) {
    # check if the column names of count_matrix match the sample names
    if (!all(colnames(count_matrix) == design_data$sample)) {
        stop("Sample names in the count matrix do not match the design data.")
    }
}

prefilter_counts <- function(count_matrix, design_data) {
    # see https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
    # getting size of smallest group
    group_sizes <- table(design_data$condition)
    smallest_group_size <- min(group_sizes)
    # keep genes with at least 10 counts over a certain number of samples
    keep <- rowSums(count_matrix >= 10) >= smallest_group_size
    filtered_count_matrix <- count_matrix[keep,]
    return(filtered_count_matrix)
}

remove_all_zero_columns <- function(df) {
    # remove columns which contains only zeros
    df <- df[, colSums(df) != 0]
    return(df)
}

replace_zero_counts_with_pseudocounts <- function(count_matrix) {
    # Add a small pseudocount of 1 to avoid zero counts
    # necessary to avoid issues with rows containing lots of (but not only) zeros
    # DESeq2 does not allow float (like 0.01) counts so we must use integers
    count_matrix[count_matrix == 0] <- 1
    return(count_matrix)
}

get_normalised_counts <- function(dds) {
    # perform normalisation
    dds <- estimateSizeFactors(dds)
    normalised_counts <- counts(dds, normalized = TRUE)
    return(normalised_counts)
}


get_cpm_counts <- function(normalised_counts, filtered_count_matrix) {
    # calculate total counts per sample (library size)
    library_sizes <- colSums(filtered_count_matrix)
    # convert normalised counts to CPM
    cpm_counts <- t(t(normalised_counts) / library_sizes * 1e6)
    # cpm_counts <- log2(cpm_counts)
    return(cpm_counts)
}

get_normalised_cpm_counts <- function(count_file, design_file) {

    print(paste('Normalizing counts in:', count_file))

    count_data <- read.csv(count_file, row.names = 1)

    # data should all be integers but sometimes they are integers converted to floats (1234 -> 1234.0)
    # DESeq2 does not accept that so we must convert them into integers
    count_data[] <- lapply(count_data, as.integer)

    design_data <- read.csv(design_file)

    count_matrix <- as.matrix(count_data)
    # in some rare datasets, columns can contain only zeros
    # we do not consider these columns
    count_matrix <- remove_all_zero_columns(count_matrix)

    # getting design data
    design_data <- design_data[design_data$sample %in% colnames(count_matrix), ]

    # check if the column names of count_matrix match the sample names
    check_samples(count_matrix, design_data)

    col_data <- data.frame(
        row.names = design_data$sample,
        condition = factor(design_data$condition)
    )

    # pre-filter genes with low counts
    filtered_count_matrix <- prefilter_counts(count_matrix, design_data)

    # if the dataframe is now empty, stop the process
    if (nrow(filtered_count_matrix) == 0) {
        message("No genes left after pre-filtering.")
        quit(save = "no", status = 100)
    }

    # add a small pseudocount to avoid zero counts
    filtered_count_matrix <- replace_zero_counts_with_pseudocounts(filtered_count_matrix)

    # create DESeq2 object
    # if the number of distinct conditions is only 1, DESeq2 returns an error
    num_unique_conditions <- length(unique(design_data$condition))
    if (num_unique_conditions == 1) {
        dds <- DESeqDataSetFromMatrix(countData = filtered_count_matrix, colData = col_data, design = ~ 1)
    } else {
        dds <- DESeqDataSetFromMatrix(countData = filtered_count_matrix, colData = col_data, design = ~ condition)
    }

    normalised_counts <- get_normalised_counts(dds)

    cpm_counts <- get_cpm_counts(normalised_counts, filtered_count_matrix)

    return(cpm_counts)
}

export_data <- function(cpm_counts, filename) {
    filename <- sub("\\.csv$", ".cpm.csv", filename)
    print(paste('Exporting normalised counts per million to:', filename))
    write.table(cpm_counts, filename, sep = ',', row.names = TRUE, col.names = NA, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

cpm_counts <- get_normalised_cpm_counts(args$count_file, args$design_file)

export_data(cpm_counts, basename(args$count_file))
