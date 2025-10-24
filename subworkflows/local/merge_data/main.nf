include { MERGE_COUNTS as MERGE_ALL_COUNTS              } from '../../../modules/local/merge/counts'
include { MERGE_COUNTS as MERGE_RNASEQ_COUNTS           } from '../../../modules/local/merge/counts'
include { MERGE_COUNTS as MERGE_MICROARRAY_COUNTS       } from '../../../modules/local/merge/counts'


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

    MERGE_RNASEQ_COUNTS ( ch_normalised_rnaseq_counts.collect() )
    MERGE_RNASEQ_COUNTS.out.counts.set { ch_merged_rnaseq_counts }

     ch_normalised_counts
        .filter { meta, file -> meta.platform == "microarray" }
        .map { meta, file -> file }
        .set { ch_normalised_microarray_counts }

    MERGE_MICROARRAY_COUNTS ( ch_normalised_microarray_counts.collect() )
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

    ch_normalised_counts
        .map {
            meta, _ -> // extracts design file and adds batch column whenever missing (for custom datasets)
                def design_content = meta.design.splitCsv( header: true )
                // if there is no batch, it is custom data
                // prepending dataset id to sample name and adding it as batch identifier
                def updated_design_content = design_content.collect { row ->
                    row.sample = row.batch ?: "custom_${meta.dataset}_${row.sample}"
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
            storeDir: "${params.outdir}/merged_datasets/"
        ) {
            item -> "${item.batch},${item.condition},${item.sample}"
        }
        .set { ch_whole_design }


    emit:
    all_counts                             = MERGE_ALL_COUNTS.out.counts
    rnaseq_counts                          = ch_merged_rnaseq_counts
    microarray_counts                      = ch_merged_microarray_counts
    whole_design                           = ch_whole_design
}
