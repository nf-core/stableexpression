include { MERGE_COUNTS as MERGE_ALL_COUNTS              } from '../../../modules/local/merge/counts'
include { MERGE_COUNTS as MERGE_RNASEQ_COUNTS           } from '../../../modules/local/merge/counts'
include { MERGE_COUNTS as MERGE_MICROARRAY_COUNTS       } from '../../../modules/local/merge/counts'
include { MERGE_DESIGNS                                 } from '../../../modules/local/merge/designs'


/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow MERGE_DATA {

    take:
    ch_normalised_counts

    main:

    // -----------------------------------------------------------------
    // MERGE COUNTS FOR EACH PLATFORM SEPARATELY
    // -----------------------------------------------------------------
    ch_normalised_counts
        .filter { meta, file -> meta.platform == "rnaseq" }
        .map { meta, file -> file }
        .set { ch_normalised_rnaseq_counts }

    MERGE_RNASEQ_COUNTS ( ch_normalised_rnaseq_counts )
    MERGE_RNASEQ_COUNTS.out.counts.set { ch_merged_rnaseq_counts }

     ch_normalised_counts
        .filter { meta, file -> meta.platform == "microarray" }
        .map { meta, file -> file }
        .set { ch_normalised_microarray_counts }

    MERGE_MICROARRAY_COUNTS ( ch_normalised_microarray_counts )
    MERGE_MICROARRAY_COUNTS.out.counts.set { ch_merged_microarray_counts }

    // -----------------------------------------------------------------
    // MERGE ALL COUNTS
    // -----------------------------------------------------------------

    ch_merged_rnaseq_counts
        .mix ( ch_merged_microarray_counts )
        .set { ch_platform_counts }

    MERGE_ALL_COUNTS( ch_platform_counts.collect())

    // -----------------------------------------------------------------
    // MERGE ALL DESIGNS IN A SINGLE TABLE
    // -----------------------------------------------------------------

    MERGE_DESIGNS(
        ch_normalised_counts.map { meta, file -> meta.design }.collect()
    )

    emit:
    all_counts                             = MERGE_ALL_COUNTS.out.counts
    rnaseq_counts                          = ch_merged_rnaseq_counts
    microarray_counts                      = ch_merged_microarray_counts
    whole_design                           = MERGE_DESIGNS.out.design
}
