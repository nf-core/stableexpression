include { GEO_GETACCESSIONS          } from '../../../modules/local/geo/getaccessions'
include { GEO_GETDATA                } from '../../../modules/local/geo/getdata'
include { addDatasetIdToMetadata     } from '../utils_nfcore_stableexpression_pipeline'
include { groupFilesByDatasetId      } from '../utils_nfcore_stableexpression_pipeline'
include { augmentMetadata            } from '../utils_nfcore_stableexpression_pipeline'
include { geoDatasetsToFetch         } from '../utils_nfcore_stableexpression_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD GEO ACCESSIONS AND DATASETS
========================================================================================
*/

workflow GEO_FETCHDATA {

    take:
    species
    skip_fetch_geo_accessions
    accessions_only
    platform
    keywords
    geo_accessions
    geo_accessions_file
    exclude_geo_accessions
    exclude_geo_accessions_file
    ch_eatlas_excluded_accessions
    ch_nb_downloaded_eatlas_datasets
    min_nb_eatlas_datasets_auto_skip_geo
    outdir


    main:

    ch_datasets = Channel.empty()
    ch_fetched_accessions = Channel.empty()

    // ------------------------------------------------------------------------------------
    // PREPARE EXCLUDED ACCESSIONS
    // ------------------------------------------------------------------------------------

    // getting accessions to exclude from GEO
    ch_eatlas_excluded_accessions
        .filter { accession -> accession.startsWith("E-GEOD-") }
        .map { accession -> accession.replace("E-GEOD-", "GSE") }
        .set { ch_excluded_eatlas_accessions }

    // parsing file listing excluded accessions
    ch_exclude_geo_accessions_file = exclude_geo_accessions_file ? Channel.fromPath(exclude_geo_accessions_file, checkIfExists: true) : Channel.empty()

    // getting accessions to exclude and preparing in the right format
    Channel.fromList( exclude_geo_accessions.tokenize(',') )
        .mix( ch_excluded_eatlas_accessions )
        .mix( ch_exclude_geo_accessions_file.splitText() )
        .unique()
        .map { acc -> acc.trim() } // removing spaces
        .set { ch_excluded_accessions }

        ch_excluded_accessions
            .collectFile(
                name: 'excluded_geo_accessions.txt',
                storeDir: "${outdir}/geo/",
                sort: true,
                newLine: true
            )
            .ifEmpty( [] )
            .set { ch_excluded_accessions_file }

    // ------------------------------------------------------------------------------------
    // GET GEO ACCESSIONS
    // ------------------------------------------------------------------------------------

    // fetching GEO accessions if applicable
    if ( !skip_fetch_geo_accessions ) {

        // checking the number of Expression Atlas datasets downloaded
        // and storing whether to skip fetching GEO accessions
        ch_geo_to_fetch = geoDatasetsToFetch( ch_nb_downloaded_eatlas_datasets, min_nb_eatlas_datasets_auto_skip_geo )

        // trick to decide whether to fetch GEO accessions or not depending on ch_geo_to_fetch
        Channel.value(species)
            .combine( ch_geo_to_fetch )
            .filter{ species_name, to_fetch -> to_fetch } // kept only when to_fetch is true
            .map { species_name, to_fetch -> species_name }
            .set { ch_species }

        // getting GEO accessions given a species name and keywords
        // keywords can be an empty string
        GEO_GETACCESSIONS(
            ch_species,
            keywords,
            platform ?: 'none',
            ch_excluded_accessions_file,
            "none"
        )

        GEO_GETACCESSIONS.out.accessions
            .splitText()
            .set { ch_fetched_accessions }

    }

    // ------------------------------------------------------------------------------------
    // PREPARE ACCESSIONS PROVIDED BY THE USER
    // ------------------------------------------------------------------------------------

    ch_geo_accessions_file = geo_accessions_file ? Channel.fromPath(geo_accessions_file, checkIfExists: true) : Channel.empty()

    Channel.fromList( geo_accessions.tokenize(',') )
        .mix( ch_geo_accessions_file.splitText() )
        .mix( ch_fetched_accessions )
        .unique()
        .filter { acc -> acc.startsWith('GSE') }
        .map { acc -> acc.trim() }
        .set { ch_accessions }

    // ------------------------------------------------------------------------------------
    // DOWNLOAD GEO DATASETS
    // ------------------------------------------------------------------------------------

    if ( !accessions_only ) {

        // Downloading GEO datasets for each accession in ch_accessions
        GEO_GETDATA(
            ch_accessions,
            species
        )

        // adding dataset id (accession + data_type) in the file meta
        // flattening in case multiple files are returned at once
        ch_design = addDatasetIdToMetadata( GEO_GETDATA.out.design.flatten() )
        ch_counts = addDatasetIdToMetadata( GEO_GETDATA.out.counts.flatten() )

        // adding design files to the meta of their respective count files
        ch_datasets = groupFilesByDatasetId( ch_design, ch_counts )

        // adding normalisation state in the meta
        augmentMetadata( ch_datasets )

    }

    emit:
    downloaded_datasets = ch_datasets

}
