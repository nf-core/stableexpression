include { FILTER_OUT_SAMPLES_WITH_TOO_MANY_ZEROS          as TOO_MANY_ZEROS                                   } from '../../../modules/local/filter_out_samples/with_too_many_zeros'
include { FILTER_OUT_SAMPLES_WITH_TOO_MANY_MISSING_VALUES as TOO_MANY_MISSING_VALUES                          } from '../../../modules/local/filter_out_samples/with_too_many_missing_values'


/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow FILTER_OUT_LOW_QUALITY_SAMPLES {

    take:
    ch_counts
    max_zero_ratio
    max_null_ratio

    main:

    // -----------------------------------------------------------------
    // REMOVE SAMPLES WITH TOO MANY ZEROS
    // -----------------------------------------------------------------

    TOO_MANY_ZEROS (
        ch_counts,
        max_zero_ratio
    )

    // -----------------------------------------------------------------
    // REMOVE SAMPLES WITH TOO MANY MISSING VALUES
    // -----------------------------------------------------------------

    TOO_MANY_MISSING_VALUES(
        TOO_MANY_ZEROS.out.counts,
        max_null_ratio
    )

    emit:
    counts                      = TOO_MANY_MISSING_VALUES.out.counts

}
