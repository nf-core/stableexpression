include { NORMALISATION_COMPUTE_CPM as COMPUTE_CPM   } from '../../../modules/local/normalisation/compute_cpm'
include { NORMALISATION_COMPUTE_TPM as COMPUTE_TPM   } from '../../../modules/local/normalisation/compute_tpm'
include { QUANTILE_NORMALISATION                     } from '../../../modules/local/quantile_normalisation'

include { GET_TRANSCRIPT_LENGTHS                     } from '../../../subworkflows/local/get_transcript_lengths'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow EXPRESSION_NORMALISATION {

    take:
    species
    ch_datasets
    normalisation_method
    quantile_norm_target_distrib
    gene_length

    main:

    //
    // MODULE: normalisation of raw count datasets (including downloaded RNA-seq datasets)
    // at the same time, removing genes that show only zero counts
    //

    ch_datasets = ch_datasets.branch {
        meta, file ->
            raw: meta.normalised == false
            normalised: meta.normalised == true
        }

    ch_raw_rnaseq_datasets_to_normalise = ch_datasets.raw.filter { meta, file -> meta.platform == 'rnaseq' }

    if ( normalisation_method == 'tpm' ) {

        if ( params.gene_length ) {

            ch_gene_length_file = channel.fromPath( params.gene_length, checkIfExists: true )

        } else {

            // download genome annotation
            // and computing length of the longest transcript gene per gene
            GET_TRANSCRIPT_LENGTHS (species)
            ch_gene_length_file = GET_TRANSCRIPT_LENGTHS.out.csv

        }

        COMPUTE_TPM(
            ch_raw_rnaseq_datasets_to_normalise,
            ch_gene_length_file
        )
        ch_raw_rnaseq_datasets_normalised = COMPUTE_TPM.out.counts

    } else { // 'cpm'

        COMPUTE_CPM( ch_raw_rnaseq_datasets_to_normalise )
        ch_raw_rnaseq_datasets_normalised = COMPUTE_CPM.out.counts

    }

    //
    // MODULE: Quantile normalisation
    //

    // putting all normalised count datasets together and performing quantile normalisation
    QUANTILE_NORMALISATION (
        ch_datasets.normalised.mix( ch_raw_rnaseq_datasets_normalised ),
        quantile_norm_target_distrib
    )


    emit:
    counts                   = QUANTILE_NORMALISATION.out.counts

}
