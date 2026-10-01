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
    # remove rows which contain only zeros
    message(paste("Nuhmber"))
    df <- df[, colSums(df, na.rm = TRUE) != 0, drop = FALSE]
    return(df)
}

remove_genes_having_only_zeros <- function(rpk) {
    # remove genes having zeros for all counts
    # it is advised to remove them before the analysis
    non_zero_rows <- rownames(rpk[apply(rpk!=0, 1, any), , drop = FALSE])
    filtered_rpk <- rpk[rownames(rpk) %in% non_zero_rows, , drop = FALSE]
    return(filtered_rpk)
}

replace_zero_counts_with_pseudocounts <- function(count_data_matrix) {
    # Add a small pseudocount of 0.01 to avoid zero counts
    count_data_matrix[count_data_matrix == 0] <- 0.01
    return(count_data_matrix)
}

get_normalised_cpm_counts <- function(count_data) {

    rpk <- as.matrix(count_data)
    # in some rare datasets, columns can contain only zeros
    # we do not consider these columns
    message("Removing columns with all zeros")
    rpk <- remove_all_zero_columns(rpk)

    if (ncol(rpk) == 0) {
        message("All columns were full of zeros.")
        write("ALL COLUMNS WERE FULL OF ZEROS", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }

    # pre-filter genes with low counts
    message("Removing samples having only zeros")
    rpk <- remove_genes_having_only_zeros(rpk)
    # if the dataframe is now empty, stop the process
    if (nrow(rpk) == 0) {
        message("No genes left after removing genes having only zeros.")
        write("NO GENES LEFT AFTER REMOVING GENES HAVING ONLY ZEROS", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }

    # Add a small pseudocount to avoid zero counts
    message("Replacing zero counts with pseudocounts")
    rpk <- replace_zero_counts_with_pseudocounts(rpk)

    message("Normalising data")
    # GeTMM: for normalization purposes, no grouping of samples
    group <- c(rep("A", ncol(rpk)))
    rpk.norm <- edgeR::DGEList(counts = rpk, group = group)

    # normalisation
    message("Calculating normalisation factors")
    rpk.norm <- edgeR::calcNormFactors(rpk.norm)

    message("Calculating CPM counts")
    norm.counts.rpk_edger <- edgeR::cpm(rpk.norm)

    message("Calculating log2(cpm +1)")
    norm.counts.rpk_edger.log2 <- log2(norm.counts.rpk_edger) + 1

    return(norm.counts.rpk_edger.log2)
}

#####################################################
# INPORT / EXPORT
#####################################################

parse_data <- function(count_file) {
    message("Parsing count file")
    count_data <- arrow::read_parquet(count_file)
    # setting gene ID as row names
    rownames(count_data) <- count_data[[1]]
    count_data <- count_data[, -1, drop = FALSE]
    return(count_data)
}

export_data <- function(count_matrix, filename) {
    filename <- sub("\\.parquet$", ".getmm.parquet", filename)
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

norm.counts.rpk_edger <- get_normalised_cpm_counts(count_data)

export_data(norm.counts.rpk_edger, basename(args$count_file))
