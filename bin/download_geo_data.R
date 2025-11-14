#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tibble"))
library(GEOquery)
library(optparse)
library(dplyr)
library(tibble)

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
        make_option("--species", type = "character", help = "Species name")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Get GEO data"
        ))
    return(args)
}


get_experiment_type <- function(data) {
  e = experimentData(data)
  experiment_type <- tolower(attr(e, "other")$type)
  if (experiment_type == "expression profiling by high throughput sequencing") {
    return("rnaseq")
  } else if (experiment_type == "expression profiling by array") {
    return("microarray")
  } else {
    return(gsub("\n", " ; ", experiment_type))
  }
}

get_platform_id <- function(data) {
  platform_id <- as.character(unique(pData(data)$platform_id))[1]
  return(platform_id)
}


format_species_name <- function(x) {
  x <- tools::toTitleCase(x)
  x <- gsub("[_-]", " ", x)
  return(x)
}


get_samples_for_species <- function(data, species) {
  pheno <- pData(data)

  # check if organism_ch2 exists
  if ("organism_ch2" %in% colnames(pheno)) {
    keep <- pheno$organism_ch1 == species & pheno$organism_ch2 == species
  } else {
    keep <- pheno$organism_ch1 == species
  }

  # return a data.frame with matching samples
  return(pheno$geo_accession[keep])
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

make_design <- function(data, accession, species) {
    metadata_df <- pData(data)
    design_df <- build_design_dataframe(metadata_df, accession)

    # get samples corresponding to species
    species_samples <- get_samples_for_species(data, species)

    # filter design dataframe
    design_df <- design_df %>%
        filter(sample %in% species_samples)

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
                quit(save = "no", status = 0)
            }
        })

    }

    return(geo_data)

}

write_warning <- function(msg) {
    file_conn <- file( WARNING_REASON_FILE, open = "a")
    cat(paste0(msg, "; "), file = file_conn, sep = "", fill = FALSE)
    close(file_conn)
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
    write_warning("RAW PROBE INTENSITIES FOUND")
  } else if (value_range[2] > 1000) {
    message("Normalized but not log-transformed (e.g. MAS5, raw intensities)")
    write_warning("PARSED INTENSITIES: NORMALIZED BUT NOT LOG-TRANSFORMED")
  } else {
    message("Unclear data origin, check GEO metadata")
    write_warning("UNCLEAR DATA ORIGIN: CHECK GEO METADATA")
  }
}


clean_count_data <- function(df) {
    message("Cleaning counts")
    # removes rows that are all NA
    df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE]
    return(df)
}


get_microarray_counts <- function(data, design_df) {
  # get count data corresponding to samples in the design
  count_df <- data.frame(exprs(data)) %>%
      select(all_of(list(design_df$sample)))
  return(count_df)
}


get_extensions <- function(file){
    extensions <- strsplit(basename(file), split="\\.")[[1]]
    return(extensions)
}


get_rnaseq_counts <- function(data, design_df) {
    pheno_data <- pData(data)

    valid_samples <- as.character(design_df$sample)
    filtered_pdata <- pheno_data[ rownames(pheno_data) %in% valid_samples, ]
    if ( nrow(filtered_pdata) == 0 ) {
        return(data.frame())
    }

    filtered_samples <- rownames(filtered_pdata)
    suppl_data <- list(filtered_pdata$supplementary_file_1)[[1]]

    count_df_list <- list()
    for (i in 1:length(filtered_samples)) {

        sample <- filtered_samples[[i]]
        url <- suppl_data[[i]]

        if ( tolower(url) == "none" || is.na(url) || url == "") {
            message(paste("Skipping sample", sample, "because no supplementary file is provided"))
            write_warning(paste("NO SUPPLEMENTARY FILE:", sample))
            next
        }

        filename <- tolower(basename(url))
        extensions <- get_extensions(filename)
        ext <- extensions[length(extensions)]
        if (ext == "gz") {
          ext <- extensions[length(extensions) - 1]
        }
        if (!(ext %in% c("txt", "tsv", "csv", "tab"))) {
          message(paste("Extension not supported:", filename))
          write_warning(paste("UNSUPPORTED EXTENSION:", ext))
          next
        }

        # skipping if it is obviously TPMs / FPKMs / RPKMs
        if (grepl("tpm", filename) | grepl("fpkm", filename) | grepl("rpkm", filename)) {
            message(paste("Skipping already normalised file", filename))
            write_warning(paste("ALREADY NORMALIZED:", filename))
            next
        }

        skip_iteration <<- FALSE
        message(paste("Downloading", filename))
        tryCatch({
            download.file(url, filename, method = "wget", quiet = TRUE)
        }, error = function(e) {
            message(paste("Unhandled error while downloading", filename, " : ", e$message))
            write_warning(paste("ERROR WHILE DOWNLOADING:", filename))
            skip_iteration <<- TRUE
        })

        # If an error occurred, skip to the next iteration
        if (skip_iteration) {
            next
        }

        separator <- NULL
        for (sep in c("\t", ",", " ")) {
            # parsing the first line to determine the separator and see if there is a header
            counts <- read.table(filename, header = FALSE, sep = sep, row.names = 1, nrows = 1)
            if (ncol(counts) > 0) {
                separator <- sep
                if (is.numeric(counts[1, 1])) {
                    has_header <- FALSE
                } else {
                    has_header <- TRUE
                }
                break
            }
        }

        if (is.null(separator)) {
            message(paste("Skipping file with no valid separator", filename))
            write_warning(paste("NO VALID SEPARATOR:", filename))
            next
        }

        counts <- read.table(filename, header = has_header, sep = separator, row.names = 1)
        # checking number of columns
        if (ncol(counts) == 1) {
            colnames(counts) <- c(sample)
        } else {
            # TODO: see how to handle multiple columns
            #colnames(counts) <- paste0(sample, "_", 1:ncol(counts))
            write_warning(paste("MULTIPLE COUNT COLUMNS:", filename))
            next
        }

        # checking type of values
        is_all_integer <- function(x) all(floor(x) == x)
        int_counts <- counts %>% select_if(is_all_integer)

        # if some values were not integers
        if (nrow(int_counts) < nrow(counts)) {
            message(paste("Skipping non-integer file", filename))
            write_warning(paste("NOT ALL INTEGERS:", filename))
            #next
        }

        # resetting the index
        counts <- tibble::rownames_to_column(counts, var = "gene_id")
        count_df_list[[i]] <- counts
    }

    # checking if all files were skipped
    if (length(count_df_list) == 0) {
        message("No valid files found")
        return(data.frame())
    }

    # full outer join
    joined_df <- Reduce(
      function(df1, df2) merge(df1, df2, by = "gene_id", all = TRUE),
      count_df_list
    )
    joined_df <- tibble::column_to_rownames(joined_df, var = "gene_id")

    return(joined_df)
}


process_data <- function(geo_data, accession, species) {

    for (i in 1:length(geo_data)) {

        data <- geo_data[[ i ]]
        file <- names(geo_data)[[ i ]]

        platform_id <- get_platform_id(data)
        experiment_type <- get_experiment_type(data)

        if ( experiment_type == "microarray") {
            message(paste("Processing microarray data:", file))

            # keeping only non empty data
            if (nrow(data) == 0) {
              write_warning(paste("NO DATA:", file))
              next
            }

            # make design dataframe
            # keep only samples corresponding to the species of interest
            design_df <- make_design(data, accession, species)

            count_df <- get_microarray_counts(data, design_df)

            # keeping only non empty data
            if (nrow(count_df) == 0 || ncol(count_df) == 0) {
              write_warning(paste("NO DATA AFTER FILTERING:", file))
              next
            }

            # checking that data are from RMA pipeline and followed proper normalisation
            check_microarray_normalisation(count_df)

        } else if (experiment_type == "rnaseq") {
            message(paste("Processing RNA-seq data:", file))
            # make design dataframe
            # keep only samples corresponding to the species of interest
            design_df <- make_design(data, accession, species)

            count_df <- get_rnaseq_counts(data, design_df)

            # keeping only non empty data
            if (nrow(count_df) == 0 || ncol(count_df) == 0) {
              message(paste("No data found for", file))
              write_warning(paste("NO DATA:", file))
              next
            }

        } else {
          message(paste("Unsupported platform:", experiment_type))
          write_warning(paste("UNSUPPORTED PLATFORM:", experiment_type))
          next
        }

        # clean counts:
        # * removes rows that are all NA
        count_df <- clean_count_data(count_df)

        # exporting count data to CSV
        export_count_data(count_df, accession, experiment_type, platform_id)

        # exporting metadata to CSV
        export_design(design_df, accession, experiment_type, platform_id)
    }
}


export_count_data <- function(count_df, accession, experiment_type, platform_id) {

    # renaming columns, to make them specific to accession and data type
    colnames(count_df) <- paste0(accession, '_', colnames(count_df))

    outfilename <- paste0(accession, '_', platform_id, '.', experiment_type, '.normalised.counts.csv')

    # exporting to CSV file
    # index represents gene names
    message(paste('Exporting count data to file', outfilename))
    write.table(count_df, outfilename, sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE)
}

export_design <- function(design_df, accession, experiment_type, platform_id) {

    new_sample_names <- paste0(accession, '_', design_df$sample)

    df <- design_df %>%
        mutate(sample = new_sample_names ) %>%
        select(sample, condition, batch)

    outfilename <- paste0(accession, '_', platform_id, '.', experiment_type, '.design.csv')
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

process_data(geo_data, args$accession, species)
