include { COMPUTE_DATASET_STATISTICS as DESCRIPTIVE_STATISTICS                     } from '../../../modules/local/compute_dataset_statistics'

/*
========================================================================================
    SUBWORKFLOW TO COMPUTE VARIOUS STATISTICS AT THE DATASET / SAMPLE LEVEL
========================================================================================
*/

workflow DATASET_ANALYSIS {

    take:
    ch_counts

    main:

    // -----------------------------------------------------------------
    // COMPUTE VARIOUS STATISTICS AT THE SAMPLE LEVEL
    // -----------------------------------------------------------------

    DESCRIPTIVE_STATISTICS ( ch_counts )


}
