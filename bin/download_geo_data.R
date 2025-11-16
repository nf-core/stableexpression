#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("stringr"))
library(GEOquery)
library(optparse)
library(dplyr)
library(tibble)
library(stringr)

options(error = traceback)

COUNT_FILE_EXTENSION <- ".counts.csv"
DESIGN_FILE_EXTENSION <- ".design.csv"
MAPPING_FILE_EXTENSION <- ".sample_name_mapping.csv"
METADATA_FILE_EXTENSION <- ".platform_metadata.csv"
BASE_REJECTED_DIR <- "rejected"

FAILURE_REASON_FILE <- "failure_reason.txt"
WARNING_REASON_FILE <- "warning_reason.txt"

#####################################################
#####################################################
# ARG PARSER
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


#####################################################
#####################################################
# UTILS
#####################################################
#####################################################

format_species_name <- function(x) {
  x <- tools::toTitleCase(x)
  x <- gsub("[_-]", " ", x)
  return(x)
}

write_warning <- function(msg) {
    message(msg)
    file_conn <- file( WARNING_REASON_FILE, open = "a")
    cat(paste0(msg, "; "), file = file_conn, sep = "", fill = FALSE)
    close(file_conn)
}


get_extensions <- function(file){
    extensions <- strsplit(basename(file), split="\\.")[[1]]
    return(extensions)
}


get_rejected_dir <- function(platform, series) {
    rejected_dir <- file.path(BASE_REJECTED_DIR, paste0(series$accession, '_', platform$id))
    dir.exists(rejected_dir) || dir.create(rejected_dir, recursive = TRUE)
    return(rejected_dir)
}

#####################################################
#####################################################
# DOWNLOAD
#####################################################
#####################################################

download_geo_data_with_retries <- function(accession, max_retries = 3, wait_time = 5) {

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

#####################################################
#####################################################
# PARSE SERIES / PLATFORM METADATA
#####################################################
#####################################################

get_experiment_data <- function(geo_data) {
    data <- geo_data[[1]]
    experiment_data <- experimentData(data)
    return(experiment_data)
}


get_experiment_type <- function(geo_data) {
  experiment_data <- get_experiment_data(geo_data)
  experiment_type <- tolower(attr(experiment_data, "other")$type)
  if (experiment_type == "expression profiling by high throughput sequencing") {
    return("rnaseq")
  } else if (experiment_type == "expression profiling by array") {
    return("microarray")
  } else {
    return(gsub("\n", " ; ", experiment_type))
  }
}


get_series_supplementary_data <- function(geo_data) {
  experiment_data <- get_experiment_data(geo_data)
  suppl_data_str <- attr(experiment_data, "other")$supplementary_file
  return(stringr::str_split(suppl_data_str, "\n")[[1]])
}


get_platform_id <- function(metadata) {
  platform_id <- as.character(unique(metadata$platform_id))[1]
  return(platform_id)
}


#####################################################
#####################################################
# RNASEQ SAMPLES
#####################################################
#####################################################

get_rnaseq_samples <- function(geo_data, design_df) {
  rnaseq_sample_df_list <- list()
  for (i in 1:length(geo_data)) {
      data <- geo_data[[ i ]]
      metadata <- pData(data)
      rnaseq_sample_df_list[[i]] <- metadata %>%
          filter(library_strategy == "RNA-Seq" & geo_accession %in% design_df$sample) %>%
          select(geo_accession)
  }
  # concatenate rows
  rnaseq_sample_df <- Reduce(
    function(df1, df2) dplyr::bind_rows(df1, df2),
    rnaseq_sample_df_list
  )
  return(rnaseq_sample_df$geo_accession)
}



#####################################################
#####################################################
# SAMPLE NAME MAPPING
#####################################################
#####################################################


make_sample_name_mapping <- function(geo_data) {
    message("Making sample name mapping")
    mapping_df_list <- list()
    for (i in 1:length(geo_data)) {
        data <- geo_data[[ i ]]
        metadata <- pData(data)
        mapping_df_list[[i]] <- metadata %>%
            mutate(
              sample_id = geo_accession,
              sample_name = title
            ) %>%
            select(sample_id, sample_name)
    }
    # concatenate rows
    mapping_df <- Reduce(
      function(df1, df2) dplyr::bind_rows(df1, df2),
      mapping_df_list
    )
    return(mapping_df)
}

rename_columns <- function(df, mapping_df) {
  id_map <- setNames(mapping_df$sample_id, mapping_df$sample_name)
  names(df) <- ifelse(
    names(df) %in% names(id_map),
    id_map[names(df)],
    names(df)
  )
  return(df)
}

#####################################################
#####################################################
# DESIGN
#####################################################
#####################################################

get_samples_for_species <- function(metadata, species) {
  # check if organism_ch2 exists
  if ("organism_ch2" %in% colnames(metadata)) {
    keep <- metadata$organism_ch1 == species & metadata$organism_ch2 == species
  } else {
    keep <- metadata$organism_ch1 == species
  }

  # return a data.frame with matching samples
  return(metadata$geo_accession[keep])
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
      mutate(sample = geo_accession) %>% # change column name geo_accession to sample
      group_by(!!!syms(columns_to_group)) %>% # group by all columns for grouping found
      mutate(group_num = cur_group_id()) %>% # create column made from group id
      ungroup() %>%
      mutate(
        condition = paste0("G", group_num), # create condition column from group number
        batch = accession
      ) %>%
      select(sample, condition, batch) %>%
      arrange(condition)

    return(design_df)
}


get_design_for_platform <- function(design_df, metadata) {
    platform_samples <- metadata$geo_accession
    platform_design_df <- design_df %>%
        filter(sample %in% platform_samples)
    return(platform_design_df)
}

get_design_for_rnaseq <- function(design_df, rnaseq_samples) {
    rnaseq_design_df <- design_df %>%
        filter(sample %in% rnaseq_samples)
    return(rnaseq_design_df)
}


make_design <- function(metadata, series) {
    design_df <- build_design_dataframe(metadata, series$accession)
    # get samples corresponding to species
    species_samples <- get_samples_for_species(metadata, series$species)
    # filter design dataframe
    design_df <- design_df %>%
        filter(sample %in% species_samples)
    return(design_df)
}


make_overall_design <- function(geo_data, series) {
    message("Making overall design")
    design_df_list <- list()
    for (i in 1:length(geo_data)) {
        data <- geo_data[[ i ]]
        metadata <- pData(data)
        # make design dataframe
        # keep only samples corresponding to the species of interest
        design_df <- make_design(metadata, series)
        design_df_list[[i]] <- design_df
    }
    # full outer join
    design_df <- Reduce(
      function(df1, df2) dplyr::bind_rows(df1, df2),
      design_df_list
    )
    return(design_df)
}


#####################################################
#####################################################
# PARSE COUNTS FROM DATA
#####################################################
#####################################################


get_microarray_counts <- function(platform) {
  # get count data corresponding to samples in the design
  counts <- data.frame(exprs(platform$data)) %>%
      select(all_of(platform$design$sample))
  return(counts)
}


get_raw_counts_from_url <- function(data_url) {

    if ( tolower(data_url) == "none" || is.na(data_url) || data_url == "") {
        write_warning(paste("MISFORMED URL:", data_url))
        return(NULL)
    }

    filename <- tolower(basename(data_url))
    extensions <- get_extensions(filename)
    ext <- extensions[length(extensions)]
    if (ext == "gz") {
      ext <- extensions[length(extensions) - 1]
    }
    if (!(ext %in% c("txt", "tsv", "csv", "tab"))) {
      write_warning(paste("UNSUPPORTED EXTENSION:", ext, "for URL:", data_url))
      return(NULL)
    }

    message(paste("Downloading", filename))
    tryCatch({
        download.file(data_url, filename, method = "wget", quiet = TRUE)
    }, error = function(e) {
        write_warning(paste("ERROR WHILE DOWNLOADING:", filename))
        return(NULL)
    })

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
        write_warning(paste("NO VALID SEPARATOR:", filename))
        return(NULL)
    }

    message(paste("Parsing", filename))
    tryCatch({
      counts <- read.table(filename, header = has_header, sep = separator, row.names = 1)
    }, error = function(e) {
        write_warning(paste("ERROR WHILE PARSING:", filename))
        return(NULL)
    })

    # removes rows that are all NA
    counts <- counts[rowSums(!is.na(counts)) > 0, , drop = FALSE]
    return(counts)
}


get_all_rnaseq_counts <- function(platform) {
    pdata <- platform$metadata
    # getting list of samples
    samples <- pdata$geo_accession
    # getting list of columns corresponding to supp data
    supplementary_cols <- grep("^supplementary_file(_\\d+)?$", names(pdata), value = TRUE)

    count_df_list <- list()
    cpt = 1
    for (i in 1:length(samples)) {
        sample <- samples[[i]]

        for (j in 1:length(supplementary_cols)) {
            data_url <- pdata[pdata$geo_accession == sample, supplementary_cols[j]]
            counts <- get_raw_counts_from_url(data_url)
            if (is.null(counts)) {
              next
            }
            # if only one column
            if (ncol(counts) == 1) {
                colnames(counts) <- c(sample)
            }
            counts <- tibble::rownames_to_column(counts, var = "gene_id")
            # adding to list
            count_df_list[[cpt]] <- counts
            cpt = cpt + 1
        }
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


#####################################################
#####################################################
# DATA QUALITY CONTROL
#####################################################
#####################################################

is_valid_microarray <- function(platform) {

  if (!all(colnames(platform$counts) %in% platform$design$sample)) {
    message("Column names do not match samples in design")
    return(FALSE)
  }

  vals <- unlist(platform$counts, use.names = FALSE)
  vals <- vals[!is.na(vals)]

  all_integers <- all(abs(vals - round(vals)) < 1e-8)
  value_range <- range(vals, na.rm = TRUE)

  if (value_range[2] <= 20) {
      message(paste(platform$id, ": normalized, log2 scale (e.g. RMA, quantile)"))
      return(TRUE)
  } else if (all_integers) {
      write_warning(paste(platform$id, ": RAW PROBE INTENSITIES FOUND"))
      return(FALSE)
  } else if (value_range[2] > 1000) {
      write_warning(paste(platform$id, ": PARSED INTENSITIES: NORMALIZED BUT NOT LOG-TRANSFORMED"))
      return(FALSE)
  } else {
    write_warning(paste(platform$id, ": UNCLEAR DATA ORIGIN: CHECK GEO METADATA"))
    return(FALSE)
  }
}

is_valid_rnaseq <- function(platform) {

  if (!all(colnames(platform$counts) %in% platform$design$sample)) {
    message(paste(platform$id, ": column names do not match samples in design"))
    return(FALSE)
  }

  # checking if all values are integers
  tryCatch({
    is_all_integer <- function(x) all(floor(x) == x)
    int_counts <- platform$counts %>% select_if(is_all_integer)
    # if some values were not integers
    if (nrow(int_counts) < nrow(platform$counts)) {
        write_warning(paste(platform$id, ": NOT ALL INTEGERS"))
        return(FALSE)
    }
  }, error = function(e) {
      write_warning(paste(platform$id, ": COULD NOT COMPUTE FLOOR"))
      return(FALSE)
  })

  return(TRUE)
}


#####################################################
#####################################################
# EXPORT
#####################################################
#####################################################

export_count_data <- function(platform, series) {
    # renaming columns, to make them specific to accession and data type
    colnames(platform$counts) <- paste0(series$accession, '_', colnames(platform$counts))

    # if nothing is left after cleaning, we still return the original data
    # so that we can have a look at it afterwards
    if (platform$type == "microarray") {
        extension <- paste0(".normalised", COUNT_FILE_EXTENSION)
    } else {
        extension <- paste0(".raw", COUNT_FILE_EXTENSION)
    }

    outfilename <- paste0(series$accession, '_', platform$id, '.', platform$type, extension)
    if (!platform$is_valid) {
        outfilename <- file.path(get_rejected_dir(platform, series), outfilename)
    }

    # exporting to CSV file
    # index represents gene names
    message(paste(platform$id, ': exporting count data to file', outfilename))
    write.table(platform$counts, outfilename, sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE)
}


export_design <- function(platform, series) {
    new_sample_names <- paste0(series$accession, '_', series$design$sample)
    design_df <- series$design %>%
        mutate(sample = new_sample_names ) %>%
        select(sample, condition, batch)

    outfilename <- paste0(series$accession, '_', platform$id, '.', platform$type, DESIGN_FILE_EXTENSION)
    if (!platform$is_valid) {
        outfilename <- file.path(get_rejected_dir(platform, series), outfilename)
    }

    message(paste(platform$id, ': exporting design data to file', outfilename))
    write.table(design_df, outfilename, sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
}


export_name_mapping <- function(platform, series) {
    outfilename <- paste0(series$accession, '_', platform$id, '.', platform$type, MAPPING_FILE_EXTENSION)
    if (!platform$is_valid) {
        outfilename <- file.path(get_rejected_dir(platform, series), outfilename)
    }
    message(paste(platform$id, ': exporting design data to file', outfilename))
    write.table(series$mapping, outfilename, sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
}

export_metadata <- function(platform, series) {
    outfilename <- paste0(series$accession, '_', platform$id, '.', platform$type, METADATA_FILE_EXTENSION)
    if (!platform$is_valid) {
        outfilename <- file.path(get_rejected_dir(platform, series), outfilename)
    }
    message(paste(platform$id, ': exporting metadata to file', outfilename))
    write.table(platform$metadata, outfilename, sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
}


#####################################################
#####################################################
# PROCESS DATA
#####################################################
#####################################################

post_process_and_export <- function(platform, series) {
    # keeping only non empty data
    if (nrow(platform$counts) == 0 || ncol(platform$counts) == 0) {
    message(paste(platform$id, ': no data found'))
      write_warning(paste(platform$id, ": NO DATA"))
      return(NULL)
    }
    # rename columns when needed
    platform$counts <- rename_columns(platform$counts, series$mapping)

    export_count_data(platform, series)
    export_design(platform, series)
    export_name_mapping(platform, series)
    export_metadata(platform, series)
}


process_platform_data <- function(platform, series) {

    platform$metadata <- pData(platform$data)
    platform$design <- get_design_for_platform(series$design, platform$metadata)
    valid_samples <- as.character(platform$design$sample)
    platform$id <- get_platform_id(platform$metadata)

    if (length(valid_samples) == 0) {
    message(paste(platform$id, ": no sample corresponding to species", series$species))
        return(NULL)
    }

    if (platform$type == "microarray") {
        platform$counts <- get_microarray_counts(platform)
        platform$is_valid <- is_valid_microarray(platform)
    } else {
        platform$counts <- get_all_rnaseq_counts(platform)
        platform$is_valid <- is_valid_rnaseq(platform)
    }

    post_process_and_export(platform, series)
}


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


main <- function() {

    args <- get_args()

    series <- list()

    series$accession <- args$accession
    series$species <- format_species_name(args$species)

    message(paste("Getting data for accession", series$accession))
    # searching and downloading expression atlas data
    geo_data <- download_geo_data_with_retries(series$accession)

    # make a single design dataframe for all samples in the series
    series$design <- make_overall_design(geo_data, series)
    if ( length(series$design) == 0 ) {
        message("No sample corresponding to species", series$species)
        write(paste("NO SAMPLES FOR SPECIES", series$species), file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    }

    # make a map associating sample names to sample IDs
    series$mapping <- make_sample_name_mapping(geo_data)

    series$experiment_type <- get_experiment_type(geo_data)

    suppl_data_urls <- get_series_supplementary_data(geo_data)
    # for now, considering suppl data as raw rnaseq data
    # TODO: check if these are always raw rnaseq data
    if (length(suppl_data_urls) > 0) {

        message("Processing supplementary data")
        for (supp_data_url in suppl_data_urls) {
            counts <- get_raw_counts_from_url(supp_data_url)
            if (is.null(counts)) {
              next
            }
            platform <- list(
                type = "rnaseq",
                id = "suppl",
                counts = counts,
                design = series$design
            )
            platform$is_valid <- is_valid_rnaseq(platform)
            post_process_and_export(platform, series)
        }

    }

    # NOTE: we consider that a series is either a microarray series OR contains RNA-seq data
    # mixed types should be found only in SuperSeries, and it is not handled for now
    if ( series$experiment_type == "microarray" ) {

        message("Processing microarray data")
        for (i in 1:length(geo_data)) {
            platform <- list(
              type = "microarray",
              data = geo_data[[ i ]]
            )
            process_platform_data(platform, series)
        }

    } else {

        rnaseq_samples <- get_rnaseq_samples(geo_data, series$design)
        if ( series$experiment_type == "rnaseq" || length(rnaseq_samples) > 0 ) {

            message("Processing RNA-seq data")
            # taking a subset of the design corresponding to bona-fide RNA-seq samples
            rnaseq_design_df <- get_design_for_rnaseq(series$design, rnaseq_samples)
            for (i in 1:length(geo_data)) {
                platform <- list(
                  type = "rnaseq",
                  data = geo_data[[ i ]]
                )
                process_platform_data(platform, series)
            }

        } else {
          write_warning(paste("UNSUPPORTED PLATFORM:", series$experiment_type))
        }
    }
}


#####################################################
# ENTRYPOINT
#####################################################
main()
