#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("dplyr"))
library(GEOquery)
library(optparse)
library(dplyr)

options(error = traceback)

FAILURE_REASON_FILE <- "failure_reason.txt"
WARNING_REASON_FILE <- "warning_reason.txt"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

get_args <- function() {
    option_list <- list(
        make_option("--accession", type = "character", help = "Accession number of GEO dataset. Example: GSE56413"),
        make_option("--species", type = "character", help = "Accession number of GEO dataset. Example: GSE56413")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Get GEO data"
        ))
    return(args)
}


format_species_name <- function(x) {
  x <- tools::toTitleCase(x)
  x <- gsub("[_-]", " ", x)
  return(x)
}


get_samples_for_species <- function(eset, species) {
  pheno <- pData(eset)

  # check if organism_ch2 exists
  if ("organism_ch2" %in% colnames(pheno)) {
    keep <- pheno$organism_ch1 == format_species_name(species) & pheno$organism_ch2 == format_species_name(species)
  } else {
    keep <- pheno$organism_ch1 == format_species_name(species)
  }

  # return a data.frame with matching samples
  pheno$geo_accession[keep]
}


get_columns_for_grouping <- function(df) {

    base_columns <- c("characteristics", "treatment_protocol", "label_protocol", "extract_protocol", "growth_protocol")

    columns_to_group <- c()
    for (base_col in base_columns) {
      ch1_col <- paste0(base_col, "_ch1")
      ch2_col <- paste0(base_col, "_ch2")

      if (ch1_col %in% colnames(df)) {
        columns_to_group <- c(columns_to_group, ch1_col)
      }
      if (ch2_col %in% colnames(df)) {
        columns_to_group <- c(columns_to_group, ch2_col)
      }
    }

    return(columns_to_group)
}


build_design_dataframe <- function(df, accession) {
    message("Build design dataframe")

    columns_to_group <- get_columns_for_grouping(df)

    design_df <- df %>%
      mutate(sample = rownames(.)) %>%
      group_by(!!!syms(columns_to_group)) %>%
      mutate(group_num = cur_group_id()) %>%
      ungroup() %>%
      mutate(
        condition = paste0("G", group_num),
        batch = accession
      ) %>%
      select(sample, condition, batch) %>%
      arrange(condition)

    return(design_df)
}


download_geo_data_with_retries <- function(accession, species, max_retries = 3, wait_time = 5) {

    success <- FALSE
    attempts <- 0

    while (!success && attempts < max_retries) {
        attempts <- attempts + 1

        tryCatch({
            geo_data <- GEOquery::getGEO( accession )
            success <- TRUE

        }, error = function(e) {

            message("Attempt ", attempts, " Message: ", e$message)

            if (attempts < max_retries) {
                warning("Retrying in ", wait_time, " seconds...")
                Sys.sleep(wait_time)

            } else {
                warning("Unhandled error: ", e$message)
                write("EXPERIMENT NOT FOUND", file = FAILURE_REASON_FILE)
            }
        })

    }

    return(geo_data)

}


check_microarray_normalisation <- function(df) {

  vals <- unlist(df, use.names = FALSE)
  vals <- vals[!is.na(vals)]

  all_integers <- all(abs(vals - round(vals)) < 1e-8)
  value_range <- range(vals, na.rm = TRUE)

  if (value_range[2] <= 20) {
    message("Normalized, log2 scale (e.g. RMA, quantile)")
  } else if (all_integers) {
    message("Raw probe intensities (unnormalized CEL-like data)")
    write("RAW PROBE INTENSITIES FOUND", file = WARNING_REASON_FILE)
  } else if (value_range[2] > 1000) {
    message("Normalized but not log-transformed (e.g. MAS5, raw intensities)")
    write("PARSED INTENSITIES: NORMALIZED BUT NOT LOG-TRANSFORMED", file = WARNING_REASON_FILE)
  } else {
    message("Unclear data origin, check GEO metadata")
    write("UNCLEAR DATA ORIGIN: CHECK GEO METADATA", file = WARNING_REASON_FILE)
  }
}


clean_count_data <- function(df) {
    message("Cleaning counts")
    # removes rows that are all NA
    df <- df[rowSums(!is.na(df)) > 0, ]

}


process_data <- function(atlas_data, accession, species) {

    eset <- geo_data[[ 1 ]]
    #print(exprs(eset))
    # Get metadata table
    metadata_df <- pData(eset)
    design_df <- build_design_dataframe(metadata_df, accession)

    # get samples corresponding to species
    species_samples <- get_samples_for_species(eset, species)

    # filter design dataframe
    design_df <- design_df %>%
        filter(sample %in% species_samples)

    if ( length(names(geo_data)) > 1 ) {
        warning("Multiple data files were found")
        write("EXPERIMENT CONTAINS MULTIPLE FILES", file = FAILURE_REASON_FILE)
    }

    file <- names(geo_data)[[ 1 ]]

    data <- geo_data [[ file ]]
    #print(fData(data))
    # get count data for samples corresponding to the species of interest
    count_df <- data.frame(exprs(data)) %>%
        select(all_of(species_samples))

    # checking that data are from RMA pipeline and followed proper normalisation
    # raises error otherwise
    check_microarray_normalisation(count_df)

    # clean counts:
    # * removes rows that are all NA
    count_df <- clean_count_data(count_df)

    # exporting count data to CSV
    export_count_data(count_df, accession)

    # exporting metadata to CSV
    export_metadata(design_df, accession)
}


export_count_data <- function(count_df, batch_id) {

    # renaming columns, to make them specific to accession and data type
    colnames(count_df) <- paste0(batch_id, '_', colnames(count_df))

    outfilename <- paste0(batch_id, '.microarray.normalised.counts.csv')

    # exporting to CSV file
    # index represents gene names
    message(paste('Exporting count data to file', outfilename))
    write.table(count_df, outfilename, sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE)
}

export_metadata <- function(design_df, batch_id) {

    new_sample_names <- paste0(batch_id, '_', design_df$sample)

    df <- design_df %>%
        mutate(sample = new_sample_names ) %>%
        select(sample, condition, batch)

    outfilename <- paste0(batch_id, '.design.csv')
    message(paste('Exporting design data to file', outfilename))
    write.table(df, outfilename, sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
}


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

cat(paste("Getting data for accession", args$accession, "\n"))

species <- format_species_name(args$species)
# searching and downloading expression atlas data
geo_data <- download_geo_data_with_retries(args$accession, species)

process_data(geo_data, args$accession, args$species)
