include { EXPRESSIONATLAS_GETDATA as EXPRESSION_ATLAS     } from '../../../modules/local/expressionatlas/getdata'
include { GEO_GETDATA as GEO                             } from '../../../modules/local/geo/getdata'

include { addDatasetIdToMetadata                        } from '../utils_nfcore_stableexpression_pipeline'
include { groupFilesByDatasetId                         } from '../utils_nfcore_stableexpression_pipeline'
include { augmentMetadata                               } from '../utils_nfcore_stableexpression_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD GEO ACCESSIONS AND DATASETS
========================================================================================
*/

workflow DOWNLOAD_PUBLIC_DATASETS {

    take:
    species
    ch_accessions


    main:

    ch_datasets = channel.empty()
    ch_fetched_accessions = channel.empty()

    ch_accessions = ch_accessions
         .branch { acc ->
             eatlas: acc.startsWith('E-')
             geo: acc.startsWith('GSE')
         }

    // ------------------------------------------------------------------------------------
    // DOWNLOAD EXPRESSION ATLAS DATASETS
    // ------------------------------------------------------------------------------------

    // Downloading Expression Atlas data for each accession in ch_accessions
    EXPRESSION_ATLAS( ch_accessions.eatlas )

    // ------------------------------------------------------------------------------------
    // DOWNLOAD GEO DATASETS
    // ------------------------------------------------------------------------------------

    // Downloading GEO datasets for each accession in ch_accessions
    GEO(
        ch_accessions.geo,
        species
    )

    ch_downloaded_counts = EXPRESSION_ATLAS.out.counts.mix ( GEO.out.counts )
    ch_downloaded_design = EXPRESSION_ATLAS.out.design.mix ( GEO.out.design )

    // adding dataset id (accession + data_type) in the file meta
    // flattening in case multiple files are returned at once
    ch_counts = addDatasetIdToMetadata( ch_downloaded_counts.flatten() )
    ch_design = addDatasetIdToMetadata( ch_downloaded_design.flatten() )

    // adding design files to the meta of their respective count files
    ch_datasets = groupFilesByDatasetId( ch_design, ch_counts )

    // adding normalisation state in the meta
    ch_datasets = augmentMetadata( ch_datasets )

    emit:
    datasets = ch_datasets

}
