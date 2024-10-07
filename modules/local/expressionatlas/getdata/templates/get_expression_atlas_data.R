#!/usr/bin/env Rscript

print(getwd())

library("ExpressionAtlas")

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

download_expression_atlas_data <- function(species, properties) {

    species <- tolower(sub("_", " ", species))

    # searching data corresponding to the provided inputs (species & properties)
    if (is.na(properties)) {
        atlasRes <- searchAtlasExperiments( species = species)
    } else {
        atlasRes <- searchAtlasExperiments( species = species, properties = properties )
    }

    # if no data were found, raising an exception
    if (length(atlasRes) == 0) {
        properties_string <- paste(properties, collapse = ", ")
        msg <- paste('No Expression Atlas experiment found with keywords', properties_string)
        stop(msg)
    }

    # if data were found, downloading them
    atlas_data <- getAtlasData( atlasRes\$Accession )
    return(atlas_data)
}

export_count_data <- function(atlas_data) {
    # looping through experiments fetched
    for (experiment_id in names(atlas_data)) {

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
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

# searching and downloading expression atlas data
atlas_data <- download_expression_atlas_data(species = '${species}', properties = '${keywords}')

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
