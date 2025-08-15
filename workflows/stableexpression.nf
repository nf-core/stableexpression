/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_FETCHDATA              } from '../subworkflows/local/expressionatlas_fetchdata'
include { IDMAPPING                              } from '../subworkflows/local/idmapping'
include { EXPRESSION_NORMALISATION               } from '../subworkflows/local/expression_normalisation'
include { MULTIQC_WORKFLOW                       } from '../subworkflows/local/multiqc'

include { MERGE_DATA                             } from '../modules/local/merge_data'
include { GENE_STATISTICS                        } from '../modules/local/gene_statistics'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STABLEEXPRESSION {

    take:
    ch_input_datasets


    main:


    ch_top_stable_genes_summary = Channel.empty()
    ch_all_genes_statistics = Channel.empty()
    ch_top_stable_genes_transposed_counts = Channel.empty()
    ch_gene_count_statistics = Channel.empty()
    ch_skewness_statistics = Channel.empty()
    ch_ks_stats = Channel.empty()
    ch_distribution_correlations = Channel.empty()

    ch_species = Channel.value( params.species.split(' ').join('_') )

    // -----------------------------------------------------------------
    // FETCH AND DOWNLOAD EXPRESSION ATLAS DATASETS IF NEEDED
    // -----------------------------------------------------------------

    EXPRESSIONATLAS_FETCHDATA( ch_species )

    if ( !params.accessions_only ) {

        // putting all datasets together (local datasets + Expression Atlas datasets)
        ch_input_datasets
            .concat( EXPRESSIONATLAS_FETCHDATA.out.downloaded_datasets )
            .set { ch_datasets }

        // -----------------------------------------------------------------
        // IDMAPPING
        // -----------------------------------------------------------------

        IDMAPPING ( ch_datasets, ch_species )

        // -----------------------------------------------------------------
        // NORMALISATION OF RAW COUNT DATASETS (INCLUDING RNA-SEQ DATASETS)
        // -----------------------------------------------------------------

        EXPRESSION_NORMALISATION(
            IDMAPPING.out.datasets,
            params.normalisation_method,
            params.quant_norm_target_distrib
        )

        EXPRESSION_NORMALISATION.out.normalised_counts.set { ch_normalised_counts }
        EXPRESSION_NORMALISATION.out.dataset_statistics.set { ch_dataset_statistics }

        // -----------------------------------------------------------------
        // MERGE COUNT FILES AND DESIGN FILES AND FILTER OUT ZERO COUNTS
        // -----------------------------------------------------------------

        MERGE_DATA(
            ch_normalised_counts.map {  meta, file -> [file] }.collect(),
            ch_dataset_statistics.map { meta, file -> [file] }.collect(),
            params.nb_top_gene_candidates
        )

        MERGE_DATA.out.candidate_gene_counts.set { ch_candidate_gene_counts }
        MERGE_DATA.out.ks_test_statistics.set { ch_ks_stats }
        MERGE_DATA.out.gene_count_statistics.set { ch_gene_count_statistics }
        MERGE_DATA.out.skewness_statistics.set { ch_skewness_statistics }
        MERGE_DATA.out.distribution_correlations.set { ch_distribution_correlations }

        // -----------------------------------------------------------------
        // GENE STATISTICS
        // -----------------------------------------------------------------

        GENE_STATISTICS(
            MERGE_DATA.out.all_counts,
            IDMAPPING.out.gene_metadata.collect(),
            IDMAPPING.out.gene_id_mapping.collect(),
            params.nb_top_gene_candidates,
            ch_ks_stats,
            params.ks_pvalue_threshold
        )

        GENE_STATISTICS.out.top_stable_genes_summary.set { ch_top_stable_genes_summary }
        GENE_STATISTICS.out.all_statistics.set { ch_all_genes_statistics }
        GENE_STATISTICS.out.top_stable_genes_transposed_counts.set { ch_top_stable_genes_transposed_counts }

    }

    // -----------------------------------------------------------------
    // MULTIQC
    // -----------------------------------------------------------------

    Channel.empty()
        .mix( ch_top_stable_genes_summary.collect() )
        .mix( ch_all_genes_statistics.collect() )
        .mix( ch_top_stable_genes_transposed_counts.collect() )
        .mix( ch_gene_count_statistics.collect() )
        .mix( ch_skewness_statistics.collect() )
        .mix( ch_ks_stats.collect() )
        .mix( ch_distribution_correlations.collect() )
        .mix( Channel.topic('all_eatlas_experiment_metadata').collect() )
        .mix( Channel.topic('filtered_eatlas_experiment_metadata').collect() )
        .set { ch_multiqc_files }

    MULTIQC_WORKFLOW( ch_multiqc_files )

    MULTIQC_WORKFLOW.out.report.toList().set { multiqc_report }


    emit:
        multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
