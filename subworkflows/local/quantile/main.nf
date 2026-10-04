include { QUANTILE_NORMALISATION as QN                     } from '../../../modules/local/normalisation/quantile'
include { MERGE_COUNTS           as MERGE_BY_PLATFORM      } from '../../../modules/local/merge_counts'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow QUANTILE {

    take:
    ch_datasets
    ch_rnaseq_whole_design
    ch_microarray_whole_design
    quantile_norm_target_distrib

    main:

    // -----------------------------------------------------------------
    // QUANTILE NORMALISATION
    // -----------------------------------------------------------------

    QN(
        ch_datasets,
        quantile_norm_target_distrib
    )
    ch_normalised_counts = QN.out.counts

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

    MERGE_BY_PLATFORM (
        ch_rnaseq.concat( ch_microarray )
    )

    ch_merged_counts_with_design = MERGE_BY_PLATFORM.out.counts
                                    .map { meta, file -> [ meta.platform, meta, file ] }
                                    .join( ch_rnaseq_whole_design.mix( ch_microarray_whole_design ) )
                                    .map { platform, meta, file, design -> [ meta, file, design ] }


    emit:
    quantiled_normalised_merged_per_platform = ch_merged_counts_with_design

}
