#!/usr/bin/env Rscript

# Written by Olivier Coen. Released under the MIT license.

suppressPackageStartupMessages(library("GEOquery"))
library(GEOquery)
library(optparse)
library(biomaRt)


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

download_geo_data_with_retries <- function(accession, species, max_retries = 3, wait_time = 5) {
    success <- FALSE
    attempts <- 0
    print(listEnsemblGenomes())
    ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")
    print(searchDatasets(ensembl_plants, pattern = species))

    while (!success && attempts < max_retries) {

        attempts <- attempts + 1
        geo_data <- GEOquery::getGEO( accession )
        eset <- geo_data[[ 1 ]]
        # inspect available sample metadata
        species_samples <- get_samples_for_species(eset, species)
        # List all variable names
        #print(colnames(pData(eset)))

        for (file in names(geo_data)) {

            data <- geo_data [[ file ]]

            #print(data)
            #counts <- exprs(data)
            #samples <- pData(data)
            #features <- fData(data)
            #print("samples")
            #print(samples)
             #print("features")
            #print(head(features))
            #print(counts)
            #print(samples)
            #print(features)

        }
        success <- TRUE

    }

    return(geo_data)
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
    cat(paste('Exporting count data to file', outfilename))
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
    cat(paste('Exporting design data to file', outfilename))
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

            if ( data_type == 'rnaseq' ) {
                result <- get_rnaseq_data(data)
            } else if ( startsWith(data_type, 'A-') ) { # typically: A-AFFY- or A-GEOD-
                result <- get_one_colour_microarray_data(data)
            } else {
                stop(paste('ERROR: Unknown data type:', data_type))
            }

        }, error = function(e) {
            print(paste("Caught an error: ", e$message))
            print(paste('ERROR: Could not get assay data for experiment ID', accession, 'and data type', data_type))
            skip_iteration <<- TRUE
        })

        # If an error occurred, skip to the next iteration
        if (skip_iteration) {
            next
        }

        #batch_id <- get_batch_id(accession, data_type)

        # exporting count data to CSV
        #export_count_data(result, batch_id)

        # exporting metadata to CSV
        #export_metadata(result, batch_id)
    }

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

# writing count data in atlas_data to specific CSV files
#process_data(atlas_data, args$accession)

