include { GEO_GETACCESSIONS          } from '../../../modules/local/geo/getaccessions'
include { GEO_GETDATA                } from '../../../modules/local/geo/getdata'
include { addDatasetIdToMetadata     } from '../utils_nfcore_stableexpression_pipeline'
include { groupFilesByDatasetId      } from '../utils_nfcore_stableexpression_pipeline'
include { augmentToMetadata          } from '../utils_nfcore_stableexpression_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD GEO ACCESSIONS AND DATASETS
========================================================================================
*/

workflow GEO_FETCHDATA {

    take:
    ch_species
    ch_excluded_accessions

    main:

    ch_datasets = Channel.empty()
    ch_fetched_accessions = Channel.empty()

    ch_geo_accessions_file = params.geo_accessions_file ? Channel.fromPath(params.geo_accessions_file, checkIfExists: true) : Channel.empty()

    Channel.fromList( params.geo_accessions.tokenize(',') )
        .mix( ch_geo_accessions_file.splitText() )
        .unique()
        .map { acc -> acc.trim() }
        .set { ch_input_accessions }

    // fetching GEO accessions if applicable
    if ( !params.skip_fetch_geo_accessions ) {

        ch_excluded_accessions
            .collectFile(
                name: 'excluded_geo_accessions.txt',
                sort: true,
                newLine: true
            )
            .ifEmpty('none')
            .set { ch_excluded_accessions_file }

        // getting GEO accessions given a species name and keywords
        // keywords can be an empty string
        def platform = params.platform?: 'none'
        GEO_GETACCESSIONS(
            ch_species,
            params.keywords,
            platform,
            ch_excluded_accessions_file,
            "none"
        )

        GEO_GETACCESSIONS.out.accessions
            .splitText()
            .set { ch_fetched_accessions }

    }

    ch_exclude_geo_accessions_file = params.exclude_geo_accessions_file ? Channel.fromPath(params.exclude_geo_accessions_file, checkIfExists: true) : Channel.empty()

    // getting accessions to exclude and preparing in the right format
    Channel.fromList( params.exclude_geo_accessions.tokenize(',') )
        .mix( ch_exclude_geo_accessions_file.splitText() )
        .unique()
        .map { acc -> acc.trim() }
        .toList()
        .map { lst -> [lst] } // list of lists : mandatory when combining in the next step
        .set { ch_excluded_accessions }

    // appending to accessions provided by the user
    // ensures that no accessions is present twice (provided by the user and fetched from GEO)
    // removing excluded accessions
    ch_input_accessions
        .mix( ch_fetched_accessions )
        .unique()
        .map { acc -> acc.trim() }
        .filter { acc -> acc.startsWith('GSE') }
        .combine ( ch_excluded_accessions )
        .filter { accession, excluded_accessions -> !(accession in excluded_accessions) }
        .map { accession, excluded_accessions -> accession }
        .set { ch_accessions }

    if ( !params.accessions_only ) {

        // Downloading GEO datasets for each accession in ch_accessions
        GEO_GETDATA(
            ch_accessions,
            ch_species
        )

        // adding dataset id (accession + data_type) in the file meta
        ch_design = addDatasetIdToMetadata( GEO_GETDATA.out.design.flatten() )
        ch_counts = addDatasetIdToMetadata( GEO_GETDATA.out.counts.flatten() )

        // adding design files to the meta of their respective count files
        ch_datasets = groupFilesByDatasetId( ch_design, ch_counts )

        // adding normalisation state in the meta
        augmentToMetadata( ch_datasets )

    }

    emit:
    downloaded_datasets = ch_datasets

}
