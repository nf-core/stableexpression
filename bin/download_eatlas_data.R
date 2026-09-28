#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

options(error = traceback)
library(optparse)
suppressPackageStartupMessages(library("SummarizedExperiment"))
library(SummarizedExperiment)

FAILURE_REASON_FILE <- "failure_reason.txt"
WARNING_REASON_FILE <- "warning_reason.txt"

EXPRESSION_ATLAS_URL_BASE <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

get_args <- function() {
    option_list <- list(
        make_option("--accession", type = "character", help = "Accession number of expression atlas experiment. Example: E-MTAB-552")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Get expression atlas data"
        ))
    return(args)
}


get_atlas_experiment <- function( accession ) {

    if( ! accession_is_valid( accession ) ) {
        stop( "Experiment accession not valid. Cannot continue." )
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Build URL
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Name of file to download
    atlasExperimentSummaryFile <- paste0(accession, "-atlasExperimentSummary.Rdata")

    # Create full URL to download R data from.
    fullUrl <- paste(EXPRESSION_ATLAS_URL_BASE, accession, atlasExperimentSummaryFile, sep = "/")

    message(paste("Downloading Expression Atlas experiment summary from:\n", fullUrl))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Download
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    connection <- url( fullUrl )

    # load into its own environment
    load_env <- new.env()

    experiment_summary <- tryCatch({

        loaded_names <- load(connection, envir = load_env)

        # there should be exactly one object in this file.
        if (length(loaded_names) != 1) {
            warning("Unexpected number of objects in Rdata file: ", length(loaded_names))
        }

        # getting the first element
        get(loaded_names[[ 1 ]], envir = load_env)

    }, error = function(e) {

        warning(e$message)
        write("ERROR OCCURED DURING DOWNLOAD", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)

    }, finally = {
        close(connection)
    })

    message(paste("Successfully downloaded experiment summary object for", accession))

    return( experiment_summary )
}

accession_is_valid <- function( accession ) {

    if( missing( accession ) ) {
        warning( "Accession missing. Cannot validate." )
        return( FALSE )
    }

    if( !grepl( "^E-\\w{4}-\\d+$", accession ) ) {
        write("EXPERIMENT ACCESSION DOES NOT LOOK LIKE AN EXPRESSION ATLAS ACCESSION", file = FAILURE_REASON_FILE)
        quit(save = "no", status = 0)
    } else {
        return( TRUE )
    }
}

download_expression_atlas_data_with_retries <- function(accession, max_retries = 3, wait_time = 5) {
    success <- FALSE
    attempts <- 0

    while (!success && attempts < max_retries) {
        attempts <- attempts + 1

        tryCatch({

            atlas_data <- get_atlas_experiment( accession )
            success <- TRUE

        }, warning = function(w) {

            # else, retrying
            message("Attempt ", attempts, " Warning: ", w$message)

            if (attempts < max_retries) {
                warning("Retrying in ", wait_time, " seconds...")
                Sys.sleep(wait_time)

            } else {

                if (grepl("550 Requested action not taken; file unavailable", w$message)) {
                    warning(w$message)
                    write("EXPERIMENT SUMMARY NOT FOUND", file = FAILURE_REASON_FILE)
                    quit(save = "no", status = 0)
                } else if (grepl("Failure when receiving data from the peer", w$message)) {
                    warning(w$message)
                    write("EXPERIMENT NOT FOUND", file = FAILURE_REASON_FILE)
                    quit(save = "no", status = 0)
                } else if (grepl("FTP status was", w$message)) {
                    warning(w$message)
                    write("FTP ERROR", file = FAILURE_REASON_FILE)
                    quit(save = "no", status = 101)
                } else {
                    warning("Unhandled warning: ", w$message)
                    write("UNKNOWN ERROR", file = FAILURE_REASON_FILE)
                    quit(save = "no", status = 0)
                }
            }

        }, error = function(e) {

            message("Attempt ", attempts, " Message: ", e$message)

            if (attempts < max_retries) {
                warning("Retrying in ", wait_time, " seconds...")
                Sys.sleep(wait_time)

            } else {

                if (grepl("Download appeared successful but no experiment summary object was found", e$message)) {
                    warning(e$message)
                    write("EXPERIMENT SUMMARY NOT FOUND", file = FAILURE_REASON_FILE)
                    quit(save = "no", status = 0)
                } else {
                    warning("Unhandled error: ", e$message)
                    write("UNKNOWN ERROR", file = FAILURE_REASON_FILE)
                    quit(save = "no", status = 0)
                }

            }
        })
    }

    return(atlas_data)
}

get_rnaseq_data <- function(data) {
    return(list(
        count_data = assays( data )$counts,
        platform = 'rnaseq',
        count_type = 'raw', # rnaseq data are raw in ExpressionAtlas
        sample_groups = colData(data)$AtlasAssayGroup
        ))
}

get_one_colour_microarray_data <- function(data) {
    return(list(
        count_data = exprs( data ),
        platform = 'microarray',
        count_type = 'normalised', # one colour microarray data are already normalised in ExpressionAtlas
        sample_groups = phenoData(data)$AtlasAssayGroup
    ))
}

get_batch_id <- function(accession, data_type) {
    batch_id <- paste0(accession, '_', data_type)
    # cleaning
    batch_id <- gsub("-", "_", batch_id)
    return(batch_id)
}

get_new_sample_names <- function(result, batch_id) {
    new_colnames <- paste0(batch_id, '_', colnames(result$count_data))
    return(new_colnames)
}

export_count_data <- function(result, batch_id) {

    # renaming columns, to make them specific to accession and data type
    colnames(result$count_data) <- get_new_sample_names(result, batch_id)

    outfilename <- paste0(batch_id, '.', result$platform, '.', result$count_type, '.counts.csv')

    # exporting to CSV file
    # index represents gene names
    message(paste('Exporting count data to file', outfilename))
    write.table(result$count_data, outfilename, sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE)
}

export_metadata <- function(result, batch_id) {

    new_colnames <- get_new_sample_names(result, batch_id)
    batch_list <- rep(batch_id, length(new_colnames))

    df <- data.frame(
        batch = batch_list,
        condition = result$sample_groups,
        sample = new_colnames
    )

    outfilename <- paste0(batch_id, '.design.csv')
    message(paste('Exporting design data to file', outfilename))
    write.table(df, outfilename, sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
}


process_data <- function(atlas_data, accession) {

    # looping through each data type (ex: 'rnaseq') in the experiment
    for (data_type in names(atlas_data)) {

        data <- atlas_data[[ data_type ]]
        skip_iteration <- FALSE

        # getting count dataframe
        tryCatch({

            if ( data_type == 'rnaseq' ) {
                result <- get_rnaseq_data(data)
            } else if ( startsWith(data_type, 'A-') ) { # typically: A-AFFY- or A-GEOD-
                result <- get_one_colour_microarray_data(data)
            } else {
                warning(paste("Unknown data type:", data_type))
                write(paste("UNKNOWN DATA TYPE:", data_type), file = WARNING_REASON_FILE, append=TRUE)
                skip_iteration <<- TRUE
            }

        }, error = function(e) {
            warning(paste("Caught an error: ", e$message))
            write(paste('ERROR: COULD NOT GET ASSAY DATA FOR EXPERIMENT ID', accession, 'AND DATA TYPE', data_type), file = WARNING_REASON_FILE, append=TRUE)
            skip_iteration <<- TRUE
        })

        # If an error occurred, skip to the next iteration
        if (skip_iteration) {
            next
        }

        batch_id <- get_batch_id(accession, data_type)

        # exporting count data to CSV
        export_count_data(result, batch_id)

        # exporting metadata to CSV
        export_metadata(result, batch_id)
    }

}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

message(paste("Getting data for accession", args$accession, "\n"))

accession <- trimws(args$accession)
if (startsWith(accession, "E-PROT")) {
    warning("Ignoring the ", accession, " experiment.")
    write("PROTEOME ACCESSIONS NOT HANDLED", file = FAILURE_REASON_FILE)
    quit(save = "no", status = 0)
}

# searching and downloading expression atlas data
atlas_data <- download_expression_atlas_data_with_retries(args$accession)

# writing count data in atlas_data to specific CSV files
process_data(atlas_data, args$accession)
