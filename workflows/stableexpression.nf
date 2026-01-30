/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GET_PUBLIC_ACCESSIONS                  } from '../subworkflows/local/get_public_accessions'
include { DOWNLOAD_PUBLIC_DATASETS               } from '../subworkflows/local/download_public_datasets'
include { ID_MAPPING                             } from '../subworkflows/local/idmapping'
include { SAMPLE_FILTERING                       } from '../subworkflows/local/sample_filtering'
include { EXPRESSION_NORMALISATION               } from '../subworkflows/local/expression_normalisation'
include { DATASET_ANALYSIS                       } from '../subworkflows/local/dataset_analysis'
include { MERGE_DATA                             } from '../subworkflows/local/merge_data'
include { GENE_STATISTICS                        } from '../subworkflows/local/gene_statistics'
include { STABILITY_SCORING                      } from '../subworkflows/local/stability_scoring'
include { MULTIQC_WORKFLOW                       } from '../subworkflows/local/multiqc'


include { AGGREGATE_RESULTS                      } from '../modules/local/aggregate_results'
include { DASH_APP                               } from '../modules/local/dash_app'

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

    ch_accessions = channel.empty()
    ch_downloaded_datasets = channel.empty()

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_most_stable_genes_summary = channel.empty()
    ch_all_genes_statistics = channel.empty()
    ch_most_stable_genes_transposed_counts = channel.empty()

    def species = params.species.split(' ').join('_').toLowerCase()

    // -----------------------------------------------------------------
    // FETCH PUBLIC ACCESSIONS
    // -----------------------------------------------------------------

    GET_PUBLIC_ACCESSIONS(
        species,
        params.skip_fetch_eatlas_accessions,
        params.fetch_geo_accessions,
        params.platform,
        params.keywords,
        channel.fromList( params.accessions.tokenize(',') ),
        params.accessions_file ? channel.fromPath(params.accessions_file, checkIfExists: true) : channel.empty(),
        channel.fromList( params.excluded_accessions.tokenize(',') ),
        params.excluded_accessions_file ? channel.fromPath(params.excluded_accessions_file, checkIfExists: true) : channel.empty(),
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

    if ( !params.accessions_only && !params.download_only ) {

        ch_counts = ch_input_datasets.mix( ch_downloaded_datasets )
        // returns an error with a message if no dataset was found
        checkCounts( ch_counts )

        // -----------------------------------------------------------------
        // IDMAPPING
        // -----------------------------------------------------------------

        // tries to map gene IDs to Ensembl IDs whenever possible
        ID_MAPPING(
            ch_counts,
            species,
            params.skip_id_mapping,
            params.skip_cleaning_gene_ids,
            params.gprofiler_target_db,
            params.gene_id_mapping,
            params.gene_metadata,
            params.min_occurrence_freq,
            params.min_occurrence_quantile,
            params.outdir
        )
        ch_counts          = ID_MAPPING.out.counts
        ch_gene_id_mapping = ID_MAPPING.out.mapping
        ch_gene_metadata   = ID_MAPPING.out.metadata
        ch_valid_gene_ids  = ID_MAPPING.out.valid_gene_ids

        // -----------------------------------------------------------------
        // FILTER OUT SAMPLES NOT VALID
        // -----------------------------------------------------------------

        SAMPLE_FILTERING (
            ch_counts,
            ch_valid_gene_ids,
            params.max_zero_ratio,
            params.max_null_ratio,
            params.outdir
        )
        ch_nb_nulls_per_sample_file = SAMPLE_FILTERING.out.nb_nulls_per_sample_file

        // -----------------------------------------------------------------
        // NORMALISATION OF RAW COUNT DATASETS (INCLUDING RNA-SEQ DATASETS)
        // -----------------------------------------------------------------

        EXPRESSION_NORMALISATION(
            species,
            SAMPLE_FILTERING.out.counts,
            params.normalisation_method,
            params.quantile_norm_target_distrib,
            params.gene_length
        )
        ch_normalised_counts = EXPRESSION_NORMALISATION.out.counts

        // -----------------------------------------------------------------
        // ANALYSIS OF NORMALISED DATASETS
        // -----------------------------------------------------------------

        DATASET_ANALYSIS(
            ch_normalised_counts
        )

        // -----------------------------------------------------------------
        // MERGE ALL DATASETS INTO ONE SINGLE DATASET
        // -----------------------------------------------------------------

        MERGE_DATA (
            ch_normalised_counts,
            ch_gene_id_mapping,
            ch_gene_metadata,
            params.missing_value_imputer,
            params.outdir
        )

        ch_all_counts      = MERGE_DATA.out.all_counts
        ch_whole_design    = MERGE_DATA.out.whole_design
        ch_platform_counts = MERGE_DATA.out.platform_counts

        // -----------------------------------------------------------------
        // COMPUTE BASE STATISTICS FOR ALL GENES
        // -----------------------------------------------------------------

        GENE_STATISTICS (
            ch_all_counts,
            ch_platform_counts,
            ch_nb_nulls_per_sample_file
        )

        ch_all_datasets_stats = GENE_STATISTICS.out.stats

        // -----------------------------------------------------------------
        // GET CANDIDATES AS REFERENCE GENE AND COMPUTES VARIOUS STABILITY VALUES
        // -----------------------------------------------------------------

        STABILITY_SCORING (
            ch_all_counts.map{ meta, file -> file },
            ch_whole_design,
            ch_all_datasets_stats,
            params.candidate_selection_descriptor,
            params.nb_top_gene_candidates,
            params.min_expr_threshold,
            params.run_genorm,
            params.stability_score_weights
        )

        ch_stats_all_genes_with_scores = STABILITY_SCORING.out.summary_statistics

        // -----------------------------------------------------------------
        // AGGREGATE ALL RESULTS FOR MULTIQC
        // -----------------------------------------------------------------

        AGGREGATE_RESULTS (
            ch_all_counts.map{ meta, file -> file }.collect(),
            ch_stats_all_genes_with_scores.collect(),
            GENE_STATISTICS.out.platform_stats.collect(),
            MERGE_DATA.out.whole_gene_metadata.collect().ifEmpty([]), // handle case where there are no mappings
            MERGE_DATA.out.whole_gene_id_mapping.collect().ifEmpty([]) // handle case where there are no mappings
        )

        ch_all_genes_summary                   = AGGREGATE_RESULTS.out.all_genes_summary
        ch_most_stable_genes_summary           = AGGREGATE_RESULTS.out.most_stable_genes_summary
        ch_most_stable_genes_transposed_counts = AGGREGATE_RESULTS.out.most_stable_genes_transposed_counts_filtered

        // -----------------------------------------------------------------
        // DASH APPLICATION
        // -----------------------------------------------------------------

        DASH_APP(
            ch_all_counts.map{ meta, file -> file }.collect(),
            ch_whole_design.collect(),
            ch_all_genes_summary.collect()
        )
        ch_versions = ch_versions.mix ( DASH_APP.out.versions )

        ch_multiqc_files = ch_multiqc_files
                            .mix( ch_most_stable_genes_summary.collect() )
                            .mix( ch_all_genes_summary.collect() )
                            .mix( ch_most_stable_genes_transposed_counts.collect() )

    }

    // -----------------------------------------------------------------
    // MULTIQC
    // -----------------------------------------------------------------

    MULTIQC_WORKFLOW(
        ch_multiqc_files,
        ch_versions,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
    )


    emit:
    multiqc_report            = MULTIQC_WORKFLOW.out.report.toList()
    most_stable_genes_summary = ch_most_stable_genes_summary

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
