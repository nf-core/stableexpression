include { MERGE_COUNTS                                      } from '../../../modules/local/merge_counts'
include { COMPUTE_GENE_STATISTICS_PER_PLATFORM              } from '../../../modules/local/compute_gene_statistics/per_platform'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow MERGE_COMPUTE_STATS_PER_PLATFORM {

    take:
    ch_normalised_counts
    platform

    main:

    // -----------------------------------------------------------------
    // MERGE COUNT FILES AND DESIGN FILES AND FILTER OUT ZERO COUNTS
    // -----------------------------------------------------------------

    MERGE_COUNTS(
        ch_normalised_counts.map  { meta, file -> [file] }.collect()
    )
    MERGE_COUNTS.out.all_counts.set { ch_merged_counts }

    // -----------------------------------------------------------------
    // GENE STATISTICS
    // -----------------------------------------------------------------

    COMPUTE_GENE_STATISTICS_PER_PLATFORM(
        ch_merged_counts,
        platform
    )

    emit:
    counts                          = ch_merged_counts
    stats                           = COMPUTE_GENE_STATISTICS_PER_PLATFORM.out.stats

}
