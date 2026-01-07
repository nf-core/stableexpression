include { MERGE_COUNTS as MERGE_PLATFORM_COUNTS         } from '../../../modules/local/merge_counts'
include { MERGE_COUNTS as MERGE_ALL_COUNTS              } from '../../../modules/local/merge_counts'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow MERGE_DATA {

    take:
    ch_normalised_counts
    ch_gene_id_mapping
    ch_gene_metadata
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

    MERGE_PLATFORM_COUNTS (
        ch_collected_rnaseq_counts.concat( ch_collected_microarray_counts )
    )

    ch_platform_counts = MERGE_PLATFORM_COUNTS.out.counts

    // -----------------------------------------------------------------
    // MERGE ALL COUNTS
    // -----------------------------------------------------------------

    ch_collected_merged_counts = ch_platform_counts
                                    .map { meta, file -> file }
                                    .collect( sort: true )
                                    .map { files -> [ [ platform: "all" ], files ] }

    MERGE_ALL_COUNTS( ch_collected_merged_counts )

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

    // -----------------------------------------------------------------
    // MERGE ALL GENE ID MAPPINGS
    // -----------------------------------------------------------------

    ch_whole_gene_id_mapping = ch_gene_id_mapping
                                .filter { it != [] } // handle case where there are no mappings
                                .splitCsv( header: true )
                                .unique()
                                .collectFile(
                                    name: 'whole_gene_id_mapping.csv',
                                    seed: "original_gene_id,gene_id",
                                    newLine: true,
                                    sort: true,
                                    storeDir: "${outdir}/idmapping/"
                                ) {
                                    item -> "${item.original_gene_id},${item.gene_id}"
                                }
                                .ifEmpty([]) // handle case where there are no mappings

    // -----------------------------------------------------------------
    // MERGE ALL GENE METADATA
    // -----------------------------------------------------------------

    ch_whole_gene_metadata = ch_gene_metadata
                                .filter { it != [] } // handle case where there are no mappings
                                .splitCsv( header: true )
                                .unique()
                                .collectFile(
                                    name: 'whole_gene_metadata.csv',
                                    seed: "gene_id,name,description",
                                    newLine: true,
                                    sort: true,
                                    storeDir: "${outdir}/idmapping/"
                                ) {
                                    item -> "${item.gene_id},${item.name},${item.description}"
                                }
                                .ifEmpty([]) // handle case where there are no mappings

    emit:
    all_counts                             = MERGE_ALL_COUNTS.out.counts
    platform_counts                        = ch_platform_counts
    whole_design                           = ch_whole_design
    whole_gene_id_mapping                  = ch_whole_gene_id_mapping
    whole_gene_metadata                    = ch_whole_gene_metadata
}
