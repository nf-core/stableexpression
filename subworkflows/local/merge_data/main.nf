include { MERGE_COUNTS as PLATFORM                      } from '../../../modules/local/merge_counts'
include { MERGE_COUNTS as GLOBAL                        } from '../../../modules/local/merge_counts'
include { IMPUTE_MISSING_VALUES                         } from '../../../modules/local/impute_missing_values'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
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


    ch_normalised_rnaseq_counts = ch_normalised_counts.filter { meta, file -> meta.platform == "rnaseq" }
    ch_normalised_microarray_counts = ch_normalised_counts.filter { meta, file -> meta.platform == "microarray" }

    ch_collected_rnaseq_counts = ch_normalised_rnaseq_counts
                                    .map { meta, file -> file }
                                    .collect( sort: true )
                                    .map { files -> [ [ platform: "rnaseq" ], files ] }

    ch_collected_microarray_counts = ch_normalised_microarray_counts
                                        .map { meta, file -> file }
                                        .collect( sort: true )
                                        .map { files -> [ [ platform: "microarray" ], files ] }

    PLATFORM (
        ch_collected_rnaseq_counts.concat( ch_collected_microarray_counts )
    )

    ch_platform_counts = PLATFORM.out.counts

    // -----------------------------------------------------------------
    // MERGE ALL COUNTS
    // -----------------------------------------------------------------

    ch_collected_merged_counts = ch_platform_counts
                                    .map { meta, file -> file }
                                    .collect( sort: true )
                                    .map { files -> [ [ platform: "all" ], files ] }

    GLOBAL( ch_collected_merged_counts.collect() )
    ch_all_counts = GLOBAL.out.counts

    // -----------------------------------------------------------------
    // IMPUTE MISSING VALUES
    // -----------------------------------------------------------------

    IMPUTE_MISSING_VALUES(
        ch_all_counts.collect(),
        missing_value_imputer
    )

    // -----------------------------------------------------------------
    // MERGE ALL DESIGNS IN A SINGLE TABLE
    // -----------------------------------------------------------------

    ch_whole_design = ch_normalised_counts
                        .map {
                            meta, file -> // extracts design file and adds batch column whenever missing (for custom datasets)
                                def design_content = meta.design.splitCsv( header: true )
                                // if there is no batch, it is custom data
                                def updated_design_content = design_content.collect { row ->
                                    row.batch = row.batch ?: "custom_${meta.dataset}"
                                    return row
                                }
                                [ updated_design_content ]
                        }
                        .flatten()
                        .unique()
                        .collectFile(
                            name: 'whole_design.csv',
                            seed: "batch,condition,sample",
                            newLine: true,
                            sort: true,
                            storeDir: "${outdir}/merged_datasets/"
                        ) {
                            item -> "${item.batch},${item.condition},${item.sample}"
                        }

    emit:
    all_imputed_counts                     = IMPUTE_MISSING_VALUES.out.counts
    all_counts                             = ch_all_counts
    platform_counts                        = ch_platform_counts
    whole_design                           = ch_whole_design
}
