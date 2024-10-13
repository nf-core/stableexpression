#!/usr/bin/env Rscript

library("ExpressionAtlas")

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

download_expression_atlas_data <- function(accession) {
    atlas_data <- getAtlasData( accession )
    return(atlas_data)
}

export_count_data <- function(atlas_data) {

    experiment_id <- names(atlas_data)[1]
    eset <- atlas_data[[ experiment_id ]]

    # looping through each data type (ex: 'rnaseq') in the experiment
    for (data_type in names(eset)) {

        data <- eset[[ data_type ]]

        skip_iteration <- FALSE
        # getting count dataframe
        tryCatch({
            df <- assays(data)\$counts
        }, error = function(e) {
            print(paste("Caught an error: ", e\$message))
            print(paste('ERROR: Could not get assay data for experiment ID', experiment_id, 'and data type', data_type))
            skip_iteration <- TRUE
        })

        # If an error occurred, skip to the next iteration
        if (skip_iteration) {
            next
        }

        outfilename <- paste0(experiment_id, '_', data_type, '.csv')
        # exporting to CSV file
        # index represents gene names
        write.table(df, outfilename, sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE)
    }

}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

# searching and downloading expression atlas data
atlas_data <- download_expression_atlas_data(accession = '${accession}')

# writing count data in atlas_data to specific CSV files
export_count_data(atlas_data)

################################################
################################################
# VERSIONS
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
expressionatlas.version <- as.character(packageVersion("ExpressionAtlas"))

version_outfile <- 'versions.yml'
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-expressionatlas:', expressionatlas.version)
    ),
    version_outfile)
