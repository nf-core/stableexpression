include { NORMALISATION_COMPUTE_CPM as COMPUTE_CPM   } from '../../../modules/local/normalisation/compute_cpm'
include { NORMALISATION_COMPUTE_CPM as COMPUTE_TPM   } from '../../../modules/local/normalisation/compute_tpm'
include { QUANTILE_NORMALISATION                     } from '../../../modules/local/quantile_normalisation'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow EXPRESSION_NORMALISATION {

    take:
    ch_datasets
    normalisation_method
    quantile_norm_target_distrib

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

    ch_datasets
        .raw.filter { meta, file -> meta.platform == 'rnaseq' }
        .map { meta, file -> [ meta, file, meta.design ] }
        .set { ch_raw_rnaseq_datasets_to_normalise }

    if ( normalisation_method == 'tpm' ) {
        COMPUTE_TPM( ch_raw_rnaseq_datasets_to_normalise )
        ch_raw_rnaseq_datasets_normalised = COMPUTE_TPM.out.counts

    } else { // 'cpm'
        COMPUTE_CPM( ch_raw_rnaseq_datasets_to_normalise )
        ch_raw_rnaseq_datasets_normalised = COMPUTE_CPM.out.counts
    }

    //
    // MODULE: Quantile normalisation
    //

    // putting all normalised count datasets together and performing quantile normalisation
    ch_datasets.normalised
        .mix( ch_raw_rnaseq_datasets_normalised )
        .set { quant_norm_input }

    QUANTILE_NORMALISATION (
        quant_norm_input,
        quantile_norm_target_distrib
    )


    emit:
    counts                   = QUANTILE_NORMALISATION.out.counts

}
