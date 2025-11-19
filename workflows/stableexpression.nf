/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_FETCHDATA              } from '../subworkflows/local/expressionatlas_fetchdata'
include { GEO_FETCHDATA                          } from '../subworkflows/local/geo_fetchdata'
include { ID_MAPPING                             } from '../subworkflows/local/idmapping'
include { EXPRESSION_NORMALISATION               } from '../subworkflows/local/expression_normalisation'
include { DATA_CLEANSING                         } from '../subworkflows/local/data_cleansing'
include { MERGE_DATA                             } from '../subworkflows/local/merge_data'
include { BASE_STATISTICS                        } from '../subworkflows/local/base_statistics'
include { STABILITY_SCORING                      } from '../subworkflows/local/stability_scoring'
include { MULTIQC_WORKFLOW                       } from '../subworkflows/local/multiqc'

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

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_top_stable_genes_summary = Channel.empty()
    ch_all_genes_statistics = Channel.empty()
    ch_top_stable_genes_transposed_counts = Channel.empty()

    def species = params.species.split(' ').join('_').toLowerCase()

    // -----------------------------------------------------------------
    // FETCH AND DOWNLOAD EXPRESSION ATLAS DATASETS IF NEEDED
    // -----------------------------------------------------------------

    EXPRESSIONATLAS_FETCHDATA( species )
    EXPRESSIONATLAS_FETCHDATA.out.downloaded_datasets.set { ch_eatlas_downloaded_datasets }

    // -----------------------------------------------------------------
    // FETCH AND DOWNLOAD GEO DATASETS IF NEEDED
    // -----------------------------------------------------------------

    GEO_FETCHDATA (
        species,
        params.skip_fetch_geo_accessions,
        params.accessions_only,
        params.platform,
        params.keywords,
        params.geo_accessions,
        params.geo_accessions_file,
        params.exclude_geo_accessions,
        params.exclude_geo_accessions_file,
        EXPRESSIONATLAS_FETCHDATA.out.accessions,
        ch_eatlas_downloaded_datasets.count(),
        params.min_nb_eatlas_datasets_auto_skip_geo,
        params.outdir
    )


    // putting all datasets together (local datasets + Expression Atlas datasets)
    ch_input_datasets
        .concat( ch_eatlas_downloaded_datasets )
        .concat( GEO_FETCHDATA.out.downloaded_datasets )
        .set { ch_counts }

    // store nb of genes and nb f samples at this stage in the meta maps
    ch_counts = storeDatasetSize( ch_counts, "nb_genes", "nb_samples" )

    // displays a message if no dataset was found
    checkCounts( ch_counts )

    if ( !params.accessions_only && !params.download_only ) {

        // -----------------------------------------------------------------
        // IDMAPPING
        // -----------------------------------------------------------------

        ch_gene_id_mapping = params.gene_id_mapping_file ? Channel.fromPath( params.gene_id_mapping_file, checkIfExists: true ) : Channel.value( [] )
        ch_gene_metadata = params.gene_metadata ? Channel.fromPath( params.gene_metadata, checkIfExists: true ) : Channel.value( [] )

        if ( !params.skip_id_mapping ) {

            // tries to map gene IDs to Ensembl IDs whenever possible
            ID_MAPPING(
                ch_counts,
                species,
                ch_gene_id_mapping,
                ch_gene_metadata
            )
            ID_MAPPING.out.counts.set { ch_counts }
            ID_MAPPING.out.mapping.set { ch_gene_id_mapping }
            ID_MAPPING.out.metadata.set { ch_gene_metadata }

        }

        ch_counts = storeDatasetSize( ch_counts, "nb_genes_after_idmapping", "nb_samples_after_idmapping" )

        // -----------------------------------------------------------------
        // NORMALISATION OF RAW COUNT DATASETS (INCLUDING RNA-SEQ DATASETS)
        // -----------------------------------------------------------------

        EXPRESSION_NORMALISATION(
            ch_counts,
            params.normalisation_method,
            params.quantile_norm_target_distrib
        )

        // -----------------------------------------------------------------
        // GET STATISTICS DATASET BY DATASET AND PERFORM SOME CLEANING OPERATIONS
        // -----------------------------------------------------------------

        DATA_CLEANSING(
            EXPRESSION_NORMALISATION.out.normalised_counts,
            params.quantile_norm_target_distrib,
            params.ks_pvalue_threshold
        )

        // -----------------------------------------------------------------
        // MERGE DATA
        // -----------------------------------------------------------------

        MERGE_DATA (
            DATA_CLEANSING.out.cleaned_counts,
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
            ch_all_counts,
            ch_stats_all_genes_with_scores,
            BASE_STATISTICS.out.rnaseq_stats.ifEmpty( [] ),
            BASE_STATISTICS.out.microarray_stats.ifEmpty( [] ),
            MERGE_DATA.out.whole_gene_metadata,
            MERGE_DATA.out.whole_gene_id_mapping
        )

        AGGREGATE_RESULTS.out.all_genes_summary.set { ch_all_genes_summary }
        AGGREGATE_RESULTS.out.top_stable_genes_summary.set { ch_top_stable_genes_summary }
        AGGREGATE_RESULTS.out.top_stable_genes_transposed_counts_filtered.set { ch_top_stable_genes_transposed_counts }

        // -----------------------------------------------------------------
        // DASH APPLICATION
        // -----------------------------------------------------------------

        DASH_APP(
            ch_all_counts,
            ch_whole_design,
            ch_all_genes_summary
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
