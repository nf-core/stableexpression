include { NORMALISATION_DESEQ2                 } from '../../../modules/local/normalisation/deseq2'
include { NORMALISATION_EDGER                  } from '../../../modules/local/normalisation/edger'
include { QUANTILE_NORMALISATION               } from '../../../modules/local/quantile_normalisation'

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

    if ( normalisation_method == 'deseq2' ) {
        NORMALISATION_DESEQ2( ch_raw_rnaseq_datasets_to_normalise )
        ch_raw_rnaseq_datasets_normalised = NORMALISATION_DESEQ2.out.cpm

    } else { // 'edger'
        NORMALISATION_EDGER( ch_raw_rnaseq_datasets_to_normalise )
        ch_raw_rnaseq_datasets_normalised = NORMALISATION_EDGER.out.cpm
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
    normalised_counts                   = QUANTILE_NORMALISATION.out.counts

}
