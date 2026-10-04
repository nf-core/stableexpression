include { MERGE_COUNTS                   } from '../../../modules/local/merge_counts'

include { mergeDesign                                   } from '../../../subworkflow/local/utils_nfcore_stableexpression_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO MERGE COUNTS AND DESIGNS
========================================================================================
*/

workflow MERGE_COUNTS_AND_DESIGNS {

    take:
    ch_datasets
    outdir
    design_filename

    main:

    // -----------------------------------------------------------------
    // MERGE ALL DESIGNS IN A SINGLE TABLE
    // -----------------------------------------------------------------

    ch_design = mergeDesign(ch_normalised_counts, outdir, design_filename)

    // -----------------------------------------------------------------
    // COLLECT ALL DATAFRAMES INTO A SINGLE SORTED LIST (FOR CONSISTENCY ACROSS RUNS)
    // -----------------------------------------------------------------

    ch_sorted_datasets = ch_datasets.map { meta, file -> file }.collect( sort: true )

    MERGE_COUNTS( ch_sorted_datasets )

    // -----------------------------------------------------------------
    // CREATING A TABLE LINKING
    // -----------------------------------------------------------------

    emit:
    counts = MERGE_COUNTS.out.counts
    design = ch_design

}
