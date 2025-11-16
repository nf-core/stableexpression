#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.
options(error = traceback)
suppressPackageStartupMessages(library("DESeq2"))
library(DESeq2)
library(optparse)

FAILURE_REASON_FILE <- "failure_reason.txt"
WARNING_REASON_FILE <- "warning_reason.txt"

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

parse_dataframe <- function(file_path, ...) {
    if (grepl("\\.csv$", file_path)) {
        data <- read.csv(file_path, ...)
    } else if (grepl("\\.tsv$", file_path)) {
        data <- read.table(file_path, sep = "\t", header = TRUE, ...)
    } else {
        write("UNSUPPORTED FILE FORMAT", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }
    return(data)
}

check_samples <- function(count_matrix, design_data) {
    # check if the column names of count_matrix match the sample names
    if (!all( colnames(count_matrix) == design_data$sample )) {
        write("SAMPLE NAMES IN COUNT MATRIX DO NOT MATCH DESIGN DATA", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }
    # check for extra samples
    extra_samples <- setdiff( colnames(count_matrix), design_data$sample )
    if (length(extra_samples) > 0) {
        write(
            "THE FOLLOWING SAMPLES ARE IN THE COUNT MATRIX BUT NOT IN DESIGN: ", paste(extra_samples, collapse = ", "),
            file = WARNING_REASON_FILE
        )
    }
}

prefilter_counts <- function(count_matrix, design_data) {
    if (ncol(count_matrix) == 1) {
        keep <- count_matrix[, 1] >= 1
    } else {
        # see https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
        # getting size of smallest group
        group_sizes <- table(design_data$condition)
        smallest_group_size <- min(group_sizes)
        # keep genes with at least 10 counts over a certain number of samples
        keep <- rowSums(count_matrix >= 1) >= smallest_group_size
    }
    filtered_count_matrix <- count_matrix[keep, , drop = FALSE] # drop = FALSE: keep dataframe structure even if only one column remains
    return(filtered_count_matrix)
}

remove_all_zero_columns <- function(df) {
    # remove columns which contains only zeros
    df <- df[, colSums(df) != 0, drop = FALSE]
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

    message("Parsing count file")
    count_data <- parse_dataframe(count_file, row.names = 1)

    # data should all be integers but sometimes they are integers converted to floats (1234 -> 1234.0)
    # DESeq2 does not accept that so we must convert them into integers
    count_data[] <- lapply(count_data, as.integer)

    count_matrix <- as.matrix(count_data)

    # in some rare datasets, columns can contain only zeros
    # we do not consider these columns
    message("Removing columns with all zeros")
    count_matrix <- remove_all_zero_columns(count_matrix)

    if (ncol(count_matrix) == 0) {
        message("All columns were full of zeros.")
        write("ALL COLUMNS WERE FULL OF ZEROS", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }

    # getting design data
    message("Parsing design file")
    design_data <- parse_dataframe(design_file)

    # removing extra samples in design table
    message("Removing extra samples in design table")
    design_data <- design_data[design_data$sample %in% colnames(count_matrix), , drop = FALSE]

    if (nrow(design_data) == 0) {
        message("Design and sample names do not match.")
        write("DESIGN AND SAMPLE NAMES DO NOT MATCH", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }

    # check if the column names of count_matrix match the sample names
    message("Checking sample names")
    check_samples(count_matrix, design_data)

    # reorder count matrix columns to match design row order
    # this is absolutely mandatory
    # see https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html at part "Count matrix input"
    count_matrix <- count_matrix[, as.character(design_data$sample), drop = FALSE]

    # pre-filter genes with low counts
    message("Pre-filtering genes")
    filtered_count_matrix <- prefilter_counts(count_matrix, design_data)

    # if the dataframe is now empty, stop the process
    if (nrow(filtered_count_matrix) == 0) {
        message("No genes left after pre-filtering.")
        write("NO GENES LEFT AFTER PRE-FILTERING", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }

    # add a small pseudocount to avoid zero counts
    message("Replacing zero counts with pseudocounts")
    filtered_count_matrix <- replace_zero_counts_with_pseudocounts(filtered_count_matrix)

    # if the number of distinct conditions is only 1, DESeq2 returns an error
    message("Creating DESeqDataSet")
    col_data <- data.frame(
        row.names = design_data$sample,
        condition = factor(design_data$condition)
    )
    num_unique_conditions <- length(unique(design_data$condition))
    if (num_unique_conditions == 1) {
        dds <- DESeqDataSetFromMatrix(countData = filtered_count_matrix, colData = col_data, design = ~ 1)
    } else {
        dds <- DESeqDataSetFromMatrix(countData = filtered_count_matrix, colData = col_data, design = ~ condition)
    }

    message("Normalising counts")
    normalised_counts <- get_normalised_counts(dds)

    message("Calculating CPM counts")
    cpm_counts <- get_cpm_counts(normalised_counts, filtered_count_matrix)

    return(cpm_counts)
}

export_data <- function(cpm_counts, filename) {
    filename <- sub("\\.(csv|tsv)$", ".cpm.csv", filename)
    message(paste('Exporting normalised counts per million to:', filename))
    write.table(cpm_counts, filename, sep = ',', row.names = TRUE, col.names = NA, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

if ( is.null(args$design_file) ) {
    message("A design dataframe must be provided.")
    quit(save = "no", status = 1)
}

message(paste("Normalising counts in", args$count_file))
cpm_counts <- get_normalised_cpm_counts(args$count_file, args$design_file)

export_data(cpm_counts, basename(args$count_file))
