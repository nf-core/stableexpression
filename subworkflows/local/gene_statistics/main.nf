include { COMPUTE_GENE_STATISTICS as GLOBAL                      } from '../../../modules/local/compute_gene_statistics'
include { COMPUTE_GENE_STATISTICS as PLATFORM                    } from '../../../modules/local/compute_gene_statistics'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow GENE_STATISTICS {

    take:
    ch_all_counts      // [ [ platform: platform, dataset_size: size], file ]
    ch_platform_counts // [ [ platform: platform, dataset_size: size], file ]
    ch_nb_nulls_per_sample_file

    main:

    // -----------------------------------------------------------------
    // PLATFORM-SPECIFIC STATISTICS
    // -----------------------------------------------------------------

    PLATFORM(
        ch_platform_counts,
        ch_nb_nulls_per_sample_file.collect()
    )


    // -----------------------------------------------------------------
    // ALL DATA
    // -----------------------------------------------------------------

    GLOBAL(
        ch_all_counts.collect(),
        ch_nb_nulls_per_sample_file.collect()
    )

    emit:
    stats                           = GLOBAL.out.stats
    platform_stats                  = PLATFORM.out.stats

}
