include { COMPUTE_GENE_STATISTICS as GLOBAL                      } from '../../../modules/local/compute_gene_statistics'
include { COMPUTE_GENE_STATISTICS as PLATFORM                    } from '../../../modules/local/compute_gene_statistics'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow GENE_STATISTICS {

    take:
    ch_all_imputed_counts
    ch_all_counts
    ch_platform_counts
    ch_ratio_nulls_per_sample_file
    max_null_ratio_valid_sample

    main:

    // -----------------------------------------------------------------
    // PLATFORM-SPECIFIC STATISTICS
    // -----------------------------------------------------------------

    // platform counts have not been imputed
    PLATFORM(
        ch_platform_counts.map{ meta, file -> [ meta, file, [] ] },
        ch_ratio_nulls_per_sample_file.collect(),
        max_null_ratio_valid_sample
    )


    // -----------------------------------------------------------------
    // ALL DATA
    // -----------------------------------------------------------------

    GLOBAL(
        ch_all_counts.join( ch_all_imputed_counts ).collect(),
        ch_ratio_nulls_per_sample_file.collect(),
        max_null_ratio_valid_sample
    )

    emit:
    stats                           = GLOBAL.out.stats
    platform_stats                  = PLATFORM.out.stats

}
