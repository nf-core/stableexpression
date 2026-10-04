include { MERGE_COUNTS                      } from '../../../modules/local/merge_counts'

include { mergeDesign                       } from '../../../subworkflow/local/utils_nfcore_stableexpression_pipeline'


/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow MICROARRAY_NORMALISATION {

    take:
    ch_datasets

    main:

    // -----------------------------------------------------------------
    // MERGE ALL DESIGNS IN A SINGLE TABLE
    // -----------------------------------------------------------------

    ch_design = mergeDesign(ch_normalised_counts, "${outdir}/merged_data/", 'microarray.original_design.csv')

    ch_sorted_datasets = ch_datasets.map { meta, file -> file }.collect( sort: true )


    MERGE_COUNTS( ch_sorted_datasets )


    ch_platform_counts = PLATFORM.out.counts




    emit:

}
