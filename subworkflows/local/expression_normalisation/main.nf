include { NORMALISATION_CPM as COMPUTE_CPM               } from '../../../modules/local/normalisation/cpm'
include { NORMALISATION_TPM as COMPUTE_TPM               } from '../../../modules/local/normalisation/tpm'
include { SCALING_NORMALISATION                          } from '../../../modules/local/scaling_normalisation'
include { QUANTILE_NORMALISATION                         } from '../../../modules/local/quantile_normalisation'

include { GET_TRANSCRIPT_LENGTHS                         } from '../../../subworkflows/local/get_transcript_lengths'
include { GETMM                                          } from '../../../subworkflows/local/getmm'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow EXPRESSION_NORMALISATION {

    take:
    species
    ch_datasets
    ch_valid_gene_ids
    skip_gene_length_normalisation
    skip_scaling
    scaling_method
    quantile_norm_target_distrib
    gff_file
    gff_url
    gene_length_file

    main:

    ch_gene_length_file = channel.empty()

    ch_datasets = ch_datasets.branch {
        meta, file ->
            rnaseq: meta.platform == 'rnaseq'
            microarray: meta.platform == 'microarray'
        }

    // ------------------------------------------------------------------------------------
    // NORMALISATION OF RNA-SEQ DATA
    // ------------------------------------------------------------------------------------

    ch_rnaseq_datasets = ch_datasets.rnaseq.branch {
        meta, file ->
            raw: meta.normalised == false
            normalised: meta.normalised == true
        }

    if ( !skip_gene_length_normalisation  ) {

        // ------------------------------------------------------------------------------------
        // GET TABLE ASSOCIATING GENE IDS TO THE LENGTH OF THEIR LARGEST TRANSCRIPT
        // ------------------------------------------------------------------------------------

        if ( gene_length_file ) {

            ch_gene_length_file = channel.fromPath( gene_length_file, checkIfExists: true )

        } else {

            // download genome annotation
            // and compute length of the longest transcript gene per gene
            GET_TRANSCRIPT_LENGTHS(
                species,
                ch_valid_gene_ids,
                gff_file,
                gff_url,
                ch_datasets.rnaseq.collect() // used to trigger this subworkflow only if RNA-seq data are present
            )
            ch_gene_length_file = GET_TRANSCRIPT_LENGTHS.out.csv

        }

        // ------------------------------------------------------------------------------------
        // NORMALISATION
        // ------------------------------------------------------------------------------------

        // normalisation on raw counts (preferred option)
        GETMM(
            ch_rnaseq_datasets.raw,
            ch_gene_length_file
        )

        // if some provided counts are already normalised, use TPM instead
        COMPUTE_TPM(
            ch_rnaseq_datasets.normalised,
            ch_gene_length_file.collect()
        )

        ch_normalised_rnaseq_datasets = GETMM.out.counts.mix( COMPUTE_TPM.out.counts )

    } else {

        // Checking that we have either raw or normalised datasets, but not both
        ch_datasets.rnaseq.collect().map { datasets ->
            def normalisation_status_list = datasets.collect{ meta, file -> meta.normalised ? 'normalised': 'raw' }.unique()
            println normalisation_status_list
        }

        COMPUTE_CPM( ch_datasets.rnaseq )
        ch_normalised_rnaseq_datasets = COMPUTE_CPM.out.counts

    }


    ch_normalised_once = ch_datasets.microarray.mix( ch_normalised_rnaseq_datasets )


    // ------------------------------------------------------------------------------------
    // SECOND NORMALISATION : SCALING
    // ------------------------------------------------------------------------------------

    //
    // put all normalised count datasets together and perform scaling (quantile normalisation or z-score)
    //

    if ( skip_scaling ) {

        ch_all_normalised = ch_normalised_once

    } else {

        if ( scaling_method == 'quantile' ) {

            QUANTILE_NORMALISATION (
                ch_normalised_once,
                quantile_norm_target_distrib
            )
            ch_all_normalised = QUANTILE_NORMALISATION.out.counts

        } else { // z-score

            SCALING_NORMALISATION( ch_normalised_once )
            ch_all_normalised = SCALING_NORMALISATION.out.counts

        }

    }

    emit:
    normalised          = ch_all_normalised
    normalised_once     = ch_normalised_once
    gene_length_file    = ch_gene_length_file

}
