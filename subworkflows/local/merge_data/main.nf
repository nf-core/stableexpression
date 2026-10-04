include { MERGE_COUNTS as PLATFORM                      } from '../../../modules/local/merge_counts'
include { MERGE_COUNTS as GLOBAL                        } from '../../../modules/local/merge_counts'
include { IMPUTE_MISSING_VALUES                         } from '../../../modules/local/impute_missing_values'

include { mergeDesign                                   } from '../../../subworkflow/local/utils_nfcore_stableexpression_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO MERGE DATASETS AND DESIGNS
========================================================================================
*/

workflow MERGE_DATA {

    take:
    ch_normalised_counts
    missing_value_imputer
    outdir

    main:

    // -----------------------------------------------------------------
    // MERGE COUNTS FOR EACH PLATFORM SEPARATELY
    // -----------------------------------------------------------------




    // -----------------------------------------------------------------
    // MERGE ALL DESIGNS IN A SINGLE TABLE
    // -----------------------------------------------------------------

    ch_whole_design = mergeDesign(ch_normalised_counts, "${outdir}/merged_data/", 'whole_design.csv')

    // -----------------------------------------------------------------
    // MERGE ALL COUNTS
    // -----------------------------------------------------------------

    ch_collected_merged_counts = ch_platform_counts
                                    .map { meta, file -> file }
                                    .collect( sort: true )
                                    .map { files -> [ [ platform: "all" ], files ] }

    GLOBAL( ch_collected_merged_counts )
    ch_all_counts = GLOBAL.out.counts

    // -----------------------------------------------------------------
    // IMPUTE MISSING VALUES
    // -----------------------------------------------------------------

    IMPUTE_MISSING_VALUES(
        ch_all_counts.collect(),
        missing_value_imputer
    )

    emit:
    all_imputed_counts                     = IMPUTE_MISSING_VALUES.out.counts
    all_counts                             = ch_all_counts
    platform_counts                        = ch_platform_counts
    whole_design                           = ch_whole_design
}
