/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GET_PUBLIC_ACCESSIONS                  } from '../subworkflows/local/get_public_accessions'
include { DOWNLOAD_PUBLIC_DATASETS               } from '../subworkflows/local/download_public_datasets'
include { ID_MAPPING                             } from '../subworkflows/local/idmapping'
include { EXPRESSION_NORMALISATION               } from '../subworkflows/local/expression_normalisation'
include { MERGE_DATA                             } from '../subworkflows/local/merge_data'
include { BASE_STATISTICS                        } from '../subworkflows/local/base_statistics'
include { STABILITY_SCORING                      } from '../subworkflows/local/stability_scoring'
include { MULTIQC_WORKFLOW                       } from '../subworkflows/local/multiqc'

include { COMPUTE_DATASET_STATISTICS             } from '../modules/local/compute_dataset_statistics'
include { AGGREGATE_RESULTS                      } from '../modules/local/aggregate_results'
include { DASH_APP                               } from '../modules/local/dash_app'

include { storeDatasetSize                       } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { checkCounts                            } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STABLEEXPRESSION {

    take:
    ch_input_datasets


    main:

    ch_accessions = Channel.empty()
    ch_downloaded_datasets = Channel.empty()

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_top_stable_genes_summary = Channel.empty()
    ch_all_genes_statistics = Channel.empty()
    ch_top_stable_genes_transposed_counts = Channel.empty()

    def species = params.species.split(' ').join('_').toLowerCase()

    // -----------------------------------------------------------------
    // FETCH PUBLIC ACCESSIONS
    // -----------------------------------------------------------------

    GET_PUBLIC_ACCESSIONS(
        species,
        params.skip_fetch_public_accessions,
        params.skip_fetch_eatlas_accessions,
        params.skip_fetch_geo_accessions,
        params.platform,
        params.keywords,
        Channel.fromList( params.accessions.tokenize(',') ),
        params.accessions_file ? Channel.fromPath(params.accessions_file, checkIfExists: true) : Channel.empty(),
        Channel.fromList( params.excluded_accessions.tokenize(',') ),
        params.excluded_accessions_file ? Channel.fromPath(params.excluded_accessions_file, checkIfExists: true) : Channel.empty(),
        params.random_sampling_size,
        params.random_sampling_seed,
        params.outdir
        )
    ch_accessions = GET_PUBLIC_ACCESSIONS.out.accessions

    // -----------------------------------------------------------------
    // DOWNLOAD GEO DATASETS IF NEEDED
    // -----------------------------------------------------------------

    if ( !params.accessions_only) {

        DOWNLOAD_PUBLIC_DATASETS (
            species,
            ch_accessions
        )
        ch_downloaded_datasets = DOWNLOAD_PUBLIC_DATASETS.out.datasets

    }

    ch_counts = ch_input_datasets.mix( ch_downloaded_datasets )

    // store nb of genes and nb f samples at this stage in the meta maps
    ch_counts = storeDatasetSize( ch_counts, "nb_genes", "nb_samples" )

    // displays a message if no dataset was found
    checkCounts( ch_counts )

    if ( !params.accessions_only && !params.download_only ) {

        // -----------------------------------------------------------------
        // IDMAPPING
        // -----------------------------------------------------------------

        // tries to map gene IDs to Ensembl IDs whenever possible
        ID_MAPPING(
            ch_counts,
            species,
            params.skip_id_mapping,
            params.gprofiler_target_db,
            params.gene_id_mapping,
            params.gene_metadata,
            params.outdir
        )
        ID_MAPPING.out.counts.set { ch_counts }
        ID_MAPPING.out.mapping.set { ch_gene_id_mapping }
        ID_MAPPING.out.metadata.set { ch_gene_metadata }

        ch_counts = storeDatasetSize( ch_counts, "nb_genes_after_idmapping", "nb_samples_after_idmapping" )

        // -----------------------------------------------------------------
        // NORMALISATION OF RAW COUNT DATASETS (INCLUDING RNA-SEQ DATASETS)
        // -----------------------------------------------------------------

        EXPRESSION_NORMALISATION(
            species,
            ch_counts,
            params.normalisation_method,
            params.quantile_norm_target_distrib
        )

        // -----------------------------------------------------------------
        // COMPUTE VARIOUS STATISTICS AT THE SAMPLE LEVEL
        // -----------------------------------------------------------------

        COMPUTE_DATASET_STATISTICS ( ch_counts )

        // -----------------------------------------------------------------
        // MERGE DATA
        // -----------------------------------------------------------------

        MERGE_DATA (
            EXPRESSION_NORMALISATION.out.counts,
            ch_gene_id_mapping,
            ch_gene_metadata
        )

        MERGE_DATA.out.all_counts.set { ch_all_counts }
        MERGE_DATA.out.whole_design.set { ch_whole_design }

        // -----------------------------------------------------------------
        // COMPUTE BASE STATISTICS FOR ALL GENES
        // -----------------------------------------------------------------

        BASE_STATISTICS (
            ch_all_counts,
            MERGE_DATA.out.rnaseq_counts,
            MERGE_DATA.out.microarray_counts
        )

        BASE_STATISTICS.out.stats.set { ch_all_datasets_stats }

        // -----------------------------------------------------------------
        // GET CANDIDATES AS REFERENCE GENE AND COMPUTES VARIOUS STABILITY VALUES
        // -----------------------------------------------------------------

        STABILITY_SCORING (
            ch_all_counts,
            ch_whole_design,
            ch_all_datasets_stats,
            params.candidate_selection_descriptor,
            params.nb_top_gene_candidates,
            params.min_expr_threshold,
            params.run_genorm,
            params.stability_score_weights
        )

        STABILITY_SCORING.out.summary_statistics.set { ch_stats_all_genes_with_scores }

        // -----------------------------------------------------------------
        // AGGREGATE ALL RESULTS FOR MULTIQC
        // -----------------------------------------------------------------

        AGGREGATE_RESULTS (
            ch_all_counts.collect(),
            ch_stats_all_genes_with_scores.collect(),
            BASE_STATISTICS.out.rnaseq_stats.ifEmpty( [] ),
            BASE_STATISTICS.out.microarray_stats.ifEmpty( [] ),
            MERGE_DATA.out.whole_gene_metadata.collect(),
            MERGE_DATA.out.whole_gene_id_mapping.collect()
        )

        AGGREGATE_RESULTS.out.all_genes_summary.set { ch_all_genes_summary }
        AGGREGATE_RESULTS.out.top_stable_genes_summary.set { ch_top_stable_genes_summary }
        AGGREGATE_RESULTS.out.top_stable_genes_transposed_counts_filtered.set { ch_top_stable_genes_transposed_counts }

        // -----------------------------------------------------------------
        // DASH APPLICATION
        // -----------------------------------------------------------------

        DASH_APP(
            ch_all_counts.collect(),
            ch_whole_design.collect(),
            ch_all_genes_summary.collect()
        )
        ch_versions = ch_versions.mix ( DASH_APP.out.versions )

        ch_multiqc_files
            .mix( ch_top_stable_genes_summary.collect() )
            .mix( ch_all_genes_summary.collect() )
            .mix( ch_top_stable_genes_transposed_counts.collect() )
            .set { ch_multiqc_files }

    }

    // -----------------------------------------------------------------
    // MULTIQC
    // -----------------------------------------------------------------

    MULTIQC_WORKFLOW(
        ch_multiqc_files,
        ch_versions
    )

    MULTIQC_WORKFLOW.out.report.toList().set { multiqc_report }


    emit:
        multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
