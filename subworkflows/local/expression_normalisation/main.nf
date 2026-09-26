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
    normalisation_method
    scaling_method
    quantile_norm_target_distrib
    gff_file
    gff_url
    gene_length_file

    main:

    // ------------------------------------------------------------------------------------
    // GET TABLE ASSOCIATING GENE IDS TO THE LENGTH OF THEIR LARGEST TRANSCRIPT
    // ------------------------------------------------------------------------------------

    ch_gene_length_file = channel.empty()

    if ( normalisation_method in ['tpm', 'getmm'] ) {

        if ( gene_length_file ) {

            ch_gene_length_file = channel.fromPath( gene_length_file, checkIfExists: true )

        } else {

            // download genome annotation
            // and computing length of the longest transcript gene per gene
            GET_TRANSCRIPT_LENGTHS(
                species,
                ch_valid_gene_ids,
                gff_file,
                gff_url
            )
            ch_gene_length_file = GET_TRANSCRIPT_LENGTHS.out.csv

        }
    }

    // ------------------------------------------------------------------------------------
    // NORMALISATION
    // ------------------------------------------------------------------------------------

    ch_datasets = ch_datasets.branch {
        meta, file ->
            raw: meta.normalised == false
            normalised: meta.normalised == true
        }

    ch_raw_rnaseq_datasets = ch_datasets.raw.filter { meta, file -> meta.platform == 'rnaseq' }


    if  ( normalisation_method == 'getmm' ) {

        GETMM(
            ch_raw_rnaseq_datasets,
            ch_gene_length_file
        )
        ch_raw_rnaseq_datasets_normalised = GETMM.out.counts

    } else if ( normalisation_method == 'tpm' ) {

        COMPUTE_TPM(
            ch_raw_rnaseq_datasets,
            ch_gene_length_file.collect()
        )
        ch_raw_rnaseq_datasets_normalised = COMPUTE_TPM.out.counts

    } else { // 'cpm'

        COMPUTE_CPM( ch_raw_rnaseq_datasets )
        ch_raw_rnaseq_datasets_normalised = COMPUTE_CPM.out.counts

    }

    ch_normalised_once = ch_datasets.normalised.mix( ch_raw_rnaseq_datasets_normalised )


    // ------------------------------------------------------------------------------------
    // SECOND NORMALISATION / SCALING
    // ------------------------------------------------------------------------------------

    //
    // put all normalised count datasets together and perform scaling (z-score) or quantile normalisation
    //

    if ( scaling_method == 'z_score' ) {

        SCALING_NORMALISATION( ch_normalised_once )
        ch_all_normalised = SCALING_NORMALISATION.out.counts

    } else { // quantile

        QUANTILE_NORMALISATION (
            ch_normalised_once,
            quantile_norm_target_distrib
        )
        ch_all_normalised = QUANTILE_NORMALISATION.out.counts

    }

    emit:
    normalised          = ch_all_normalised
    normalised_once     = ch_normalised_once
    gene_length_file    = ch_gene_length_file

}
