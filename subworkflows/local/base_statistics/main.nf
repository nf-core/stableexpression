include { COMPUTE_BASE_STATISTICS as COMPUTE_GLOBAL_STATISTICS             } from '../../../modules/local/compute_base_statistics'
include { COMPUTE_BASE_STATISTICS as COMPUTE_PLATFORM_STATISTICS           } from '../../../modules/local/compute_base_statistics'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow BASE_STATISTICS {

    take:
    ch_all_counts      // [ [ platform: platform, dataset_size: size], file ]
    ch_platform_counts // [ [ platform: platform, dataset_size: size], file ]

    main:

    // -----------------------------------------------------------------
    // PLATFORM-SPECIFIC STATISTICS
    // -----------------------------------------------------------------

    COMPUTE_PLATFORM_STATISTICS( ch_platform_counts )


    // -----------------------------------------------------------------
    // ALL DATA
    // -----------------------------------------------------------------

    COMPUTE_GLOBAL_STATISTICS( ch_all_counts )

    emit:
    stats                           = COMPUTE_GLOBAL_STATISTICS.out.stats
    platform_stats                  = COMPUTE_PLATFORM_STATISTICS.out.stats

}
