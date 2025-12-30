include { REMOVE_SAMPLES_NOT_VALID                          } from '../../../modules/local/remove_samples_not_valid'


/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow FILTER_DATASETS {

    take:
    ch_counts

    main:

    // -----------------------------------------------------------------
    // REMOVE SAMPLES WITH TOO MANY ZEROS
    // -----------------------------------------------------------------

    REMOVE_SAMPLES_NOT_VALID ( ch_counts )

    emit:
    counts                      = REMOVE_SAMPLES_NOT_VALID.out.counts

}
