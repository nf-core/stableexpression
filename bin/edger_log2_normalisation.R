#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

library(edgeR)
library(optparse)
library(arrow)

FAILURE_REASON_FILE <- "failure_reason.txt"
WARNING_REASON_FILE <- "warning_reason.txt"

#####################################################
# ARGPARSER
#####################################################

get_args <- function() {

    option_list <- list(
        optparse::make_option("--counts", dest = 'count_file', help = "Path to input count file")
    )

    args <- optparse::parse_args(optparse::OptionParser(
        option_list = option_list,
        description = "Normalize counts using edgeR"
        ))

    return(args)
}

#####################################################
# COMPUTE NORMALISATION
#####################################################

remove_all_zero_columns <- function(df) {
    # remove samples (columns) having only zeros
    message("Removing samples with all zeros")
    df[, colSums(df, na.rm = TRUE) != 0, drop = FALSE]
}

remove_genes_having_only_zeros <- function(rpk) {
    # remove genes (rows) having only zeros
    rpk[rowSums(rpk != 0, na.rm = TRUE) > 0, , drop = FALSE]
}

get_normalised_cpm_counts <- function(count_data) {

    counts <- as.matrix(count_data)

    message("Removing samples with all zeros")
    counts <- remove_all_zero_columns(counts)
    if (ncol(counts) == 0) {
        message("All columns were full of zeros.")
        write("ALL COLUMNS WERE FULL OF ZEROS", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }

    message("Removing genes having only zeros")
    counts <- remove_genes_having_only_zeros(counts)
    if (nrow(counts) == 0) {
        message("No genes left after removing genes having only zeros.")
        write("NO GENES LEFT AFTER REMOVING GENES HAVING ONLY ZEROS", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }

    message("Normalising data")
    # GeTMM: for normalization purposes, no grouping of samples
    dge <- edgeR::DGEList(counts = counts, group = rep("A", ncol(counts)))

    message("Calculating normalisation factors")
    dge <- edgeR::calcNormFactors(dge)

    message("Calculating log2(cpm + 1)")
    log2(edgeR::cpm(dge) + 1)
}

#####################################################
# INPORT / EXPORT
#####################################################

parse_data <- function(count_file) {
    message("Parsing count file")
    count_data <- as.data.frame(arrow::read_parquet(count_file))
    # setting gene ID as row names
    rownames(count_data) <- count_data[[1]]
    count_data[, -1, drop = FALSE]
}

export_data <- function(count_matrix, filename) {
    filename <- sub("\\.parquet$", ".edger_log2.parquet", filename)
    message(paste('Exporting normalised data to:', filename))
    # putting row names (gene ids) back in one column
    df <- data.frame(
      gene_id = rownames(count_matrix),
      count_matrix,
      row.names = NULL,
      check.names = FALSE
    )
    arrow::write_parquet(df, filename)
}

#####################################################
# MAIN
#####################################################

args <- get_args()

message(paste('Normalizing counts in:', args$count_file))

count_data <- parse_data(args$count_file)

normalised_counts <- get_normalised_cpm_counts(count_data)

export_data(normalised_counts, basename(args$count_file))
