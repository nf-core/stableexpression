#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

suppressPackageStartupMessages(library("affy"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("dplyr"))

# Load library
library(affy)
library(optparse)
library(AnnotationDbi)
library(dplyr)
library(tibble)

options(error = traceback)

# we need to install the affy package manually while disabling threading
# when installed through conda, we get: ERROR; return code from pthread_create() is 22
if (!requireNamespace("affy", quietly = TRUE)) {
    BiocManager::install("affy", configure.args="--disable-threading", force = TRUE, quiet = TRUE)
}


#####################################################
#####################################################
# ARG PARSER
#####################################################
#####################################################

get_args <- function() {
    option_list <- list(
        make_option("--input", help = "Folder containing CEL files"),
        make_option("--target-gene-id-db", dest = "target_gene_id_db", help = "Target database for gene IDs (ENSEMBL or ENTREZID)")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Normalize microarray data using RMA"
    ))
    return(args)
}

get_probe_id_mapping <- function(data, annot_db, target_gene_id_db, stringent) {

    probe_ids <- rownames(data)
    annotations <- AnnotationDbi::select(
        annot_db,
        keys = probe_ids,
        columns = c(target_gene_id_db),
        keytype = "PROBEID"
    )

    if (stringent) {
        annotations <- annotations %>%
            group_by(PROBEID) %>%
            filter(n_distinct(.data[[target_gene_id_db]], na.rm = TRUE) == 1) %>%
            ungroup()
    }

    return(annotations)
}

replace_probe_ids_by_target_ids <- function(data, annotations, target_gene_id_db) {
    data <- as.data.frame(data)
    data$PROBEID <- rownames(data)

    data <- merge(annotations, data, by = "PROBEID", all.x = TRUE)

    # computing mean of probe values for each gene
    data <- data %>%
        group_by(.data[[target_gene_id_db]]) %>%
        summarise(across(where(is.numeric), function(x) mean(x, na.rm = TRUE))) %>%
        ungroup()

    data <- tibble::column_to_rownames(data, var = target_gene_id_db)
    return(data)
}


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


main <- function() {

    args <- get_args()

    # Read CEL files from a directory
    message("Reading CEL files from", args$input)
    data <- ReadAffy(celfile.path = args$input)

    message("Installing annotation database")
    db_name <- paste0(annotation(data), ".db")
    if (!requireNamespace(db_name, quietly = TRUE)) {
        BiocManager::install(db_name, quiet = TRUE)
    }
    library(db_name, character.only = TRUE)

    # Normalize using RMA (most common method)
    eset <- rma(data)
    # Extract normalized expression values
    message("Extracting normalized expression values")
    normalised_data <- exprs(eset)

    annotations <- get_probe_id_mapping(
      normalised_data,
      annot_db = get(db_name), # Get the database object using get()
      target_gene_id_db = args$target_gene_id_db,
      stringent = TRUE
    )

    normalised_data_df <- replace_probe_ids_by_target_ids(normalised_data, annotations, args$target_gene_id_db)

    # cleaning colnames
    colnames(normalised_data_df) <- sub("\\..*", "", colnames(normalised_data_df))

    # Save results
    message("Saving results to normalised_expression.csv")
    write.csv(normalised_data, "normalised_expression.csv")

}

main()
