#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

# Get variation coefficient from count data for each gene
# The variation coefficient is the ratio of the standard deviation to the mean
# We want genes that are neither expressed too much nor too little
# Metadata (name and description) are used to annotate the results
# Likewise, mappings (original gene ids) are used to better associate gene ids with their original gene ids
# Usage:
# get_variation_coefficients.R --counts <count_files> --metadata <metadata_files> --mappings <mapping_files>


library(optparse)
library(data.table)
library(dplyr)
library(tibble)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


get_args <- function() {

    option_list <- list(
        make_option("--counts", dest = 'count_files', help = "Count files to join"),
        make_option("--metadata", dest = 'metadata_files', help = "Metadata files to concatenate"),
        make_option("--mappings", dest = 'mapping_files', help = "Mapping files to concatenate")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Get variation coefficient from count data for each gene"
        ))

    return(args)
}

merge_count_files <- function(file_list) {
    # Read CSV files
    dfs <- lapply(file_list, read.csv, row.names = 1, header=TRUE, stringsAsFactors = FALSE)
    # Join dataframes
    concat_df <- dfs[[1]]
    for (df in dfs[-1]) {
        # Outer join on row names (includes all rows from both data frames)
        concat_df <- merge(concat_df, df, by = "row.names", all = TRUE)
        # Reset row names
        rownames(concat_df) <- concat_df$Row.names
        concat_df <- concat_df[, -1]
    }
    return(concat_df)
}


handle_na_values <- function(df) {
    # Replace NA values with 0 (the genes were not present so it is equivalent to a 0 count)
    df[is.na(df)] <- 0
    return(df)
}


concat_and_remove_duplicates <- function(file_list) {
    # Read CSV files
    dfs <- lapply(file_list, read.csv, row.names = NULL, header = TRUE, stringsAsFactors = FALSE)
    combined_df <- rbindlist(dfs, use.names = TRUE)
    # Remove duplicates
    combined_df <- distinct(combined_df)
    return(combined_df)
}

aggregate_mappings <- function(mapping_df) {
    # group by new gene IDs and gets the list of distinct original gene IDs for each group
    aggregated_df <- mapping_df %>%
        group_by(new) %>%
        summarise(original_gene_ids = list(unique(original)), .groups = "drop")
    # convert the list column to a string representation
    # separate the original gene IDs with a semicolon
    aggregated_df$original_gene_ids <- sapply(aggregated_df$original_gene_ids, function(x) paste(x, collapse = ";"))
    return(aggregated_df)
}


average_log2 <- function(row) {
    # the dataframe has already been filtered to exclude rows where mean is 0
    return(mean(log2(row + 1))) # adds 1 to avoid log(0) and to stabilize variance
}

get_variation_coefficient <- function(count_data) {

    print('Getting coefficients of variation')

    # handling NA values (genes that are not found in all datasets)
    count_data <- handle_na_values(count_data)

    # we want genes that are neither expressed too much nor too little
    # filter the dataframe to exclude rows where row mean is in the top 5% or bottom 5%
    # determine the percentile thresholds
    row_means <- rowMeans(count_data)
    lower_threshold <- quantile(row_means, 0.05)
    upper_threshold <- quantile(row_means, 0.95)

    count_data <- count_data[row_means >= lower_threshold & row_means <= upper_threshold, ]

    # calculate the coefficient of variation (cv)
    # as the ratio of the standard deviation to the mean
    row_means <- rowMeans(count_data)
    row_sds <- apply(count_data, 1, sd)
    cv <- row_sds / row_means

    av_log_cpm <- apply(count_data, 1, average_log2)
    # combine results into a dataframe
    df <- data.frame(variation_coefficient = cv, average_log_cpm = av_log_cpm)
    # order dataframe (from lowest to highest variation coefficient)
    df <- df[order(df$variation_coefficient, decreasing = FALSE), ]

    return(df)
}

merge_data <- function(cv_df, metadata_df, mapping_df) {
    # Convert rownames to a column for joining
    cv_df <- cv_df %>% rownames_to_column("index")

    # merge with metadata
    cv_df <- cv_df %>%
        left_join(metadata_df, by = c("index" = "gene_id"))
    # merge with mappings
    cv_df <- cv_df %>%
        left_join(mapping_df, by = c("index" = "new"))

    # convert back to data frame with rownames if needed
    cv_df <- cv_df %>% column_to_rownames("index")
    return(cv_df)
}

export_data <- function(cv_df, count_data) {

    count_outfilename <- 'all_normalised_counts.csv'
    print(paste('Exporting normalised counts to:', count_outfilename))
    write.table(count_data, sep=",", file=count_outfilename, row.names = TRUE, col.names = NA, quote = FALSE)

    cv_outfilename <- 'variation_coefficients.csv'
    print(paste('Exporting variation coefficients to:', cv_outfilename))
    write.table(cv_df, sep=",", file=cv_outfilename, row.names = TRUE, col.names = NA, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

# getting variation coefficient for each gene
count_file_list <- strsplit(args$count_files, " ")[[1]]
count_data <- merge_count_files(count_file_list)
cv_df <- get_variation_coefficient(count_data)

# associating gene ids with metadata (name and description)
metadata_file_list <- strsplit(args$metadata_files, " ")[[1]]
metadata_df <- concat_and_remove_duplicates(metadata_file_list)

# associating gene ids with their original gene ids
mapping_file_list <- strsplit(args$mapping_files, " ")[[1]]
mapping_df <- concat_and_remove_duplicates(mapping_file_list)
mapping_df <- aggregate_mappings(mapping_df)

# merging gene ids / coefficients of variation, etc. with name, description and original gene ids
cv_df <- merge_data(cv_df, metadata_df, mapping_df)

export_data(cv_df, count_data)
