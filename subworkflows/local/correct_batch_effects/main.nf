include { MERGE_COUNTS as MERGE_BY_PLATFORM      } from '../../../modules/local/merge_counts'
include { RECOMBAT                               } from '../../../modules/local/recombat'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow CORRECT_BATCH_EFFECTS {

    take:
    ch_datasets
    ch_rnaseq_whole_design
    ch_microarray_whole_design

    main:

    // -----------------------------------------------------------------
    // MERGE ALL DATASETS TOGETHER FOR EACH PLATFORM SEPARATELY
    // -----------------------------------------------------------------

    ch_rnaseq_datasets     = ch_normalised_counts.filter { meta, file -> meta.platform == "rnaseq" }
    ch_microarray_datasets = ch_normalised_counts.filter { meta, file -> meta.platform == "microarray" }

    ch_rnaseq = ch_rnaseq_datasets
                    .map { meta, file -> file }
                    .collect( sort: true )
                    .map { files -> [ [ platform: "rnaseq" ], files ] }

    ch_microarray = ch_microarray_datasets
                        .map { meta, file -> file }
                        .collect( sort: true )
                        .map { files -> [ [ platform: "microarray" ], files ] }

    MERGE_BY_PLATFORM( ch_rnaseq.concat( ch_microarray ) )

    ch_counts_merged_by_platform = MERGE_BY_PLATFORM.out.counts

    // -----------------------------------------------------------------
    // RECOMBAT
    // -----------------------------------------------------------------

    ch_merged_counts_with_design = ch_counts_merged_by_platform
                                    .map { meta, file -> [ meta.platform, meta, file ] }
                                    .join( ch_rnaseq_whole_design.mix( ch_microarray_whole_design ) )
                                    .map { platform, meta, file, design -> [ meta, file, design ] }

    RECOMBAT( ch_merged_counts_with_design )


    emit:
    corrected_per_platform = RECOMBAT.out.counts


}
