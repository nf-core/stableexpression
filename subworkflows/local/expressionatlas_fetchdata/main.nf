include { EXPRESSIONATLAS_GETACCESSIONS          } from '../../../modules/local/expressionatlas/getaccessions'
include { EXPRESSIONATLAS_GETDATA                } from '../../../modules/local/expressionatlas/getdata'
include { addDatasetIdToMetadata     } from '../utils_nfcore_stableexpression_pipeline'
include { groupFilesByDatasetId      } from '../utils_nfcore_stableexpression_pipeline'
include { augmentToMetadata          } from '../utils_nfcore_stableexpression_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow EXPRESSIONATLAS_FETCHDATA {

    take:
    ch_species


    main:

    ch_eatlas_datasets = Channel.empty()
    ch_fetched_accessions = Channel.empty()

    ch_eatlas_accessions_file = params.eatlas_accessions_file ? Channel.fromPath(params.eatlas_accessions_file, checkIfExists: true) : Channel.empty()

    Channel.fromList( params.eatlas_accessions.tokenize(',') )
        .mix( ch_eatlas_accessions_file.splitText() )
        .unique()
        .map { it -> it.trim() }
        .set { ch_input_accessions }

    // fetching Expression Atlas accessions if applicable
    if ( !params.skip_fetch_eatlas_accessions || params.keywords ) {

        // getting Expression Atlas accessions given a species name and keywords
        // keywords can be an empty string
        def platform = params.platform?: 'none'
        EXPRESSIONATLAS_GETACCESSIONS(
            ch_species,
            params.keywords,
            platform
        )

        EXPRESSIONATLAS_GETACCESSIONS.out.accessions
            .splitText()
            .set { ch_fetched_accessions }

    }

    ch_exclude_eatlas_accessions_file = params.exclude_eatlas_accessions_file ? Channel.fromPath(params.exclude_eatlas_accessions_file, checkIfExists: true) : Channel.empty()

    // getting accessions to exclude and preparing in the right format
    Channel.fromList( params.exclude_eatlas_accessions.tokenize(',') )
        .mix( ch_exclude_eatlas_accessions_file.splitText() )
        .unique()
        .map { it -> it.trim() }
        .toList()
        .map { lst -> [lst] } // list of lists : mandatory when combining in the next step
        .set { ch_excluded_accessions }

    // appending to accessions provided by the user
    // ensures that no accessions is present twice (provided by the user and fetched from E. Atlas)
    // removing E-PROT- accessions
    // removing excluded accessions
    ch_input_accessions
        .mix( ch_fetched_accessions )
        .unique()
        .map { it -> it.trim() }
        .filter { it.startsWith('E-') && !it.startsWith('E-PROT-') }
        .combine ( ch_excluded_accessions )
        .filter { accession, excluded_accessions -> !(accession in excluded_accessions) }
        .map { accession, excluded_accessions -> accession }
        .set { ch_accessions }

    if ( !params.accessions_only ) {

        // Downloading Expression Atlas data for each accession in ch_accessions
        EXPRESSIONATLAS_GETDATA( ch_accessions )

        // adding dataset id (accession + data_type) in the file meta
        ch_design = addDatasetIdToMetadata( EXPRESSIONATLAS_GETDATA.out.design.flatten() )
        ch_counts = addDatasetIdToMetadata( EXPRESSIONATLAS_GETDATA.out.counts.flatten() )

        // adding design files to the meta of their respective count files
        ch_eatlas_datasets = groupFilesByDatasetId( ch_design, ch_counts )

        // adding normalisation state in the meta
        augmentToMetadata( ch_eatlas_datasets )

    }

    emit:
    downloaded_datasets = ch_eatlas_datasets
    accessions          = ch_accessions

}
