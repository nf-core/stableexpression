include { COMPUTE_BASE_STATISTICS                                                   } from '../../../modules/local/compute_base_statistics'
include { COMPUTE_BASE_STATISTICS as COMPUTE_BASE_STATISTICS_FOR_RNASEQ             } from '../../../modules/local/compute_base_statistics'
include { COMPUTE_BASE_STATISTICS as COMPUTE_BASE_STATISTICS_FOR_MICROARRAY         } from '../../../modules/local/compute_base_statistics'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow BASE_STATISTICS {

    take:
    ch_all_counts
    ch_rnaseq_counts
    ch_microarray_counts

    main:

    // -----------------------------------------------------------------
    // PLATFORM-SPECIFIC STATISTICS
    // -----------------------------------------------------------------

    COMPUTE_BASE_STATISTICS_FOR_RNASEQ(
        ch_rnaseq_counts,
        "rnaseq"
    )

    COMPUTE_BASE_STATISTICS_FOR_MICROARRAY(
        ch_microarray_counts,
        "microarray"
    )

    // -----------------------------------------------------------------
    // ALL DATA
    // -----------------------------------------------------------------

    COMPUTE_BASE_STATISTICS (
        ch_all_counts,
        []
    )

    emit:
    stats                           = COMPUTE_BASE_STATISTICS.out.stats
    rnaseq_stats                    = COMPUTE_BASE_STATISTICS_FOR_RNASEQ.out.stats
    microarray_stats                = COMPUTE_BASE_STATISTICS_FOR_MICROARRAY.out.stats

}
