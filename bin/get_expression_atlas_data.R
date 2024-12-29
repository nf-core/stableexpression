#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

library(ExpressionAtlas)
library(optparse)

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

download_expression_atlas_data_with_retries <- function(accession, max_retries = 3, wait_time = 5) {
    success <- FALSE
    attempts <- 0

    while (!success && attempts < max_retries) {
    attempts <- attempts + 1
    tryCatch({
        atlas_data <- getAtlasData( accession )
        success <- TRUE
        message("Download successful on attempt ", attempts)
    }, error = function(e) {
        message("Attempt ", attempts, " failed: ", e$message)
        if (attempts < max_retries) {
        message("Retrying in ", wait_time, " seconds...")
        Sys.sleep(wait_time)
        } else {
        message("All attempts failed. Please check the URL or your connection.")
        }
    })
    }

    return(atlas_data)
}

get_rnaseq_data <- function(data) {

    return(list(
        count_data = assays( data )$counts,
        count_type = 'raw',
        sample_groups = colData(data)$AtlasAssayGroup
        ))
}

get_one_colour_microarray_data <- function(data) {

    return(list(
        count_data = exprs( data ),
        count_type = 'normalized',
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

    outfilename <- paste0(batch_id, '.', result$count_type, '.csv')

    # exporting to CSV file
    # index represents gene names
    print(paste('Exporting count data to file', outfilename))
    write.table(result$count_data, outfilename, sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE)
}

export_metadata <- function(result, batch_id) {

    new_colnames <- get_new_sample_names(result, batch_id)
    batch_list <- rep(batch_id, length(new_colnames))

    df <- data.frame(batch = batch_list, condition = result$sample_groups, sample = new_colnames)

    outfilename <- paste0(batch_id, '.design.csv')
    print(paste('Exporting metadata to file', outfilename))
    write.table(df, outfilename, sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)
}


process_data <- function(atlas_data, accession) {

    eset <- atlas_data[[ accession ]]

    # looping through each data type (ex: 'rnaseq') in the experiment
    for (data_type in names(eset)) {

        data <- eset[[ data_type ]]

        skip_iteration <- FALSE
        # getting count dataframe
        tryCatch({

            if (data_type == 'rnaseq') {
                result <- get_rnaseq_data(data)
            } else if (startsWith(data_type, 'A-AFFY-')) {
                result <- get_one_colour_microarray_data(data)
            } else {
                stop(paste('ERROR: Unknown data type:', data_type))
            }

        }, error = function(e) {
            print(paste("Caught an error: ", e$message))
            print(paste('ERROR: Could not get assay data for experiment ID', accession, 'and data type', data_type))
            skip_iteration <- TRUE
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

# searching and downloading expression atlas data
atlas_data <- download_expression_atlas_data_with_retries(args$accession)

# writing count data in atlas_data to specific CSV files
process_data(atlas_data, args$accession)

