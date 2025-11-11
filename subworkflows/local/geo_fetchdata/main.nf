include { GEO_GETACCESSIONS          } from '../../../modules/local/geo/getaccessions'
include { GEO_GETDATA                } from '../../../modules/local/geo/getdata'
include { addDatasetIdToMetadata     } from '../utils_nfcore_stableexpression_pipeline'
include { groupFilesByDatasetId      } from '../utils_nfcore_stableexpression_pipeline'
include { augmentMetadata          } from '../utils_nfcore_stableexpression_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD GEO ACCESSIONS AND DATASETS
========================================================================================
*/

workflow GEO_FETCHDATA {

    take:
    species
    ch_eatlas_excluded_accessions

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
    ch_exclude_geo_accessions_file = params.exclude_geo_accessions_file ? Channel.fromPath(params.exclude_geo_accessions_file, checkIfExists: true) : Channel.empty()

    // getting accessions to exclude and preparing in the right format
    Channel.fromList( params.exclude_geo_accessions.tokenize(',') )
        .mix( ch_excluded_eatlas_accessions )
        .mix( ch_exclude_geo_accessions_file.splitText() )
        .unique()
        .map { acc -> acc.trim() } // removing spaces
        .set { ch_excluded_accessions }

        ch_excluded_accessions
            .collectFile(
                name: 'excluded_geo_accessions.txt',
                storeDir: "${params.outdir}/geo/",
                sort: true,
                newLine: true
            )
            .ifEmpty('none')
            .set { ch_excluded_accessions_file }

    // ------------------------------------------------------------------------------------
    // GET GEO ACCESSIONS
    // ------------------------------------------------------------------------------------

    // fetching GEO accessions if applicable
    if ( !params.skip_fetch_geo_accessions ) {

        // getting GEO accessions given a species name and keywords
        // keywords can be an empty string
        def platform = params.platform ?: 'none'
        GEO_GETACCESSIONS(
            species,
            params.keywords,
            platform,
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

    ch_geo_accessions_file = params.geo_accessions_file ? Channel.fromPath(params.geo_accessions_file, checkIfExists: true) : Channel.empty()

    Channel.fromList( params.geo_accessions.tokenize(',') )
        .mix( ch_geo_accessions_file.splitText() )
        .mix( ch_fetched_accessions )
        .unique()
        .filter { acc -> acc.startsWith('GSE') }
        .map { acc -> acc.trim() }
        .set { ch_accessions }

    // ------------------------------------------------------------------------------------
    // DOWNLOAD GEO DATASETS
    // ------------------------------------------------------------------------------------

    if ( !params.accessions_only ) {

        // Downloading GEO datasets for each accession in ch_accessions
        GEO_GETDATA(
            ch_accessions,
            species
        )

        // adding dataset id (accession + data_type) in the file meta
        ch_design = addDatasetIdToMetadata( GEO_GETDATA.out.design )
        ch_counts = addDatasetIdToMetadata( GEO_GETDATA.out.counts )

        // adding design files to the meta of their respective count files
        ch_datasets = groupFilesByDatasetId( ch_design, ch_counts )

        // adding normalisation state in the meta
        augmentMetadata( ch_datasets )

    }

    emit:
    downloaded_datasets = ch_datasets

}
