include { MERGE_COUNTS                           } from '../../../modules/local/merge/counts'
include { MERGE_DESIGNS                          } from '../../../modules/local/merge/designs'
include { COMPUTE_GLOBAL_GENE_STATISTICS         } from '../../../modules/local/compute_gene_statistics/global'

include { MERGE_COMPUTE_STATS_PER_PLATFORM as MERGE_COMPUTE_STATS_MICROARRAY } from '../merge_compute_stats_per_platform'
include { MERGE_COMPUTE_STATS_PER_PLATFORM as MERGE_COMPUTE_STATS_RNASEQ     } from '../merge_compute_stats_per_platform'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow MERGE_COMPUTE_STATS {

    take:
    ch_normalised_counts
    ch_gene_metadata
    ch_gene_id_mapping

    main:

    MERGE_COMPUTE_STATS_RNASEQ (
        ch_normalised_counts.filter { meta, file -> meta.platform == "rnaseq" },
        "rnaseq"
    )

    MERGE_COMPUTE_STATS_MICROARRAY (
        ch_normalised_counts.filter { meta, file -> meta.platform == "microarray" },
        "microarray"
    )

    // -----------------------------------------------------------------
    // MERGE RNASEQ AND MICROARRAY COUNTS
    // -----------------------------------------------------------------

    Channel.empty()
        .mix ( MERGE_COMPUTE_STATS_RNASEQ.out.counts )
        .mix ( MERGE_COMPUTE_STATS_MICROARRAY.out.counts )
        .set { ch_all_counts }

    MERGE_COUNTS( ch_all_counts.collect() )

    // -----------------------------------------------------------------
    // MERGE ALL DESIGNS IN A SINGLE TABLE
    // -----------------------------------------------------------------
    ch_normalised_counts
        .map { meta, file -> meta.design }
        .set { ch_designs }

    MERGE_DESIGNS( ch_designs.collect() )

    // -----------------------------------------------------------------
    // GENE STATISTICS
    // -----------------------------------------------------------------

    Channel.empty()
        .mix ( MERGE_COMPUTE_STATS_RNASEQ.out.stats )
        .mix ( MERGE_COMPUTE_STATS_MICROARRAY.out.stats )
        .set { ch_platform_statistics }

    COMPUTE_GLOBAL_GENE_STATISTICS(
        MERGE_COUNTS.out.all_counts,
        ch_platform_statistics.collect(),
        ch_gene_metadata.collect(),
        ch_gene_id_mapping.collect(),
        params.nb_top_gene_candidates
    )

    emit:
    top_stable_genes_summary               = COMPUTE_GLOBAL_GENE_STATISTICS.out.top_stable_genes_summary
    all_genes_statistics                   = COMPUTE_GLOBAL_GENE_STATISTICS.out.all_statistics
    top_stable_genes_transposed_counts     = COMPUTE_GLOBAL_GENE_STATISTICS.out.top_stable_genes_transposed_counts
}
