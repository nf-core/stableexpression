include { DATASET_STATISTICS                   } from '../../../modules/local/dataset_statistics'
include { CLEAN_COUNT_DATA                     } from '../../../modules/local/clean_count_data'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow DATA_CLEANSING {

    take:
    ch_quantile_normalised_datasets
    quantile_normalisation_target_distribution
    ks_pvalue_threshold

    main:

    //
    // Get global stats for each sample in each dataset
    //

    DATASET_STATISTICS(
        ch_quantile_normalised_datasets,
        quantile_normalisation_target_distribution
    )

    //
    // Filter out aberrant samples and perform some sorting / cleaning
    //

    CLEAN_COUNT_DATA (
        ch_quantile_normalised_datasets.join( DATASET_STATISTICS.out.stats ),
        ks_pvalue_threshold
    )


    emit:
    cleaned_counts                   = CLEAN_COUNT_DATA.out.counts

}
