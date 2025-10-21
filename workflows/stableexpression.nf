/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_FETCHDATA              } from '../subworkflows/local/expressionatlas_fetchdata'
include { GEO_FETCHDATA                          } from '../subworkflows/local/geo_fetchdata'
include { IDMAPPING                              } from '../subworkflows/local/idmapping'
include { EXPRESSION_NORMALISATION               } from '../subworkflows/local/expression_normalisation'
include { DATA_CLEANSING                         } from '../subworkflows/local/data_cleansing'
include { MERGE_DATA                             } from '../subworkflows/local/merge_data'
include { BASE_STATISTICS                        } from '../subworkflows/local/base_statistics'
include { STABILITY_SCORING                      } from '../subworkflows/local/stability_scoring'
include { MULTIQC_WORKFLOW                       } from '../subworkflows/local/multiqc'

include { AGGREGATE_RESULTS                      } from '../modules/local/aggregate_results'
include { DASH_APP                               } from '../modules/local/dash_app'

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

    ch_species = Channel.value( params.species.split(' ').join('_') )

    // -----------------------------------------------------------------
    // FETCH AND DOWNLOAD EXPRESSION ATLAS DATASETS IF NEEDED
    // -----------------------------------------------------------------

    EXPRESSIONATLAS_FETCHDATA( ch_species )

    // getting accessions to exclude from GEO
    EXPRESSIONATLAS_FETCHDATA.out.accessions
        .filter { accession -> accession.startsWith("E-GEOD-") }
        .map { accession -> accession.replace("E-GEOD-", "GSE")}
        .set { ch_excluded_geo_accessions }

    GEO_FETCHDATA (
        ch_species,
        ch_excluded_geo_accessions
    )

    if ( !params.accessions_only && !params.download_only ) {

        // putting all datasets together (local datasets + Expression Atlas datasets)
        ch_input_datasets
            .concat( EXPRESSIONATLAS_FETCHDATA.out.downloaded_datasets )
            .concat( GEO_FETCHDATA.out.downloaded_datasets )
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
            params.quantile_normalisation_target_distribution
        )

        // -----------------------------------------------------------------
        // GET STATISTICS DATASET BY DATASET AND PERFORM SOME CLEANING OPERATIONS
        // -----------------------------------------------------------------

        DATA_CLEANSING(
            EXPRESSION_NORMALISATION.out.normalised_counts,
            params.quantile_normalisation_target_distribution,
            params.ks_pvalue_threshold
        )

        // -----------------------------------------------------------------
        // MERGE DATA
        // -----------------------------------------------------------------

        MERGE_DATA ( DATA_CLEANSING.out.cleaned_counts )
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

        // -----------------------------------------------------------------
        // GET CANDIDATES AS REFERENCE GENE AND COMPUTES VARIOUS STABILITY VALUES
        // -----------------------------------------------------------------

        STABILITY_SCORING (
            ch_all_counts,
            ch_whole_design,
            BASE_STATISTICS.out.stats
        )

        // -----------------------------------------------------------------
        // AGGREGATE ALL RESULTS FOR MULTIQC
        // -----------------------------------------------------------------

        AGGREGATE_RESULTS (
            ch_all_counts,
            STABILITY_SCORING.out.summary_statistics,
            IDMAPPING.out.gene_metadata,
            IDMAPPING.out.gene_id_mapping
        )

        AGGREGATE_RESULTS.out.top_stable_genes_summary.set { ch_top_stable_genes_summary }
        AGGREGATE_RESULTS.out.stats_all_genes.set { ch_all_genes_statistics }
        AGGREGATE_RESULTS.out.top_stable_genes_transposed_counts_filtered.set { ch_top_stable_genes_transposed_counts }

    }

    // -----------------------------------------------------------------
    // MULTIQC
    // -----------------------------------------------------------------

    Channel.empty()
        .mix( ch_top_stable_genes_summary.collect() )
        .mix( ch_all_genes_statistics.collect() )
        .mix( ch_top_stable_genes_transposed_counts.collect() )
        .mix( Channel.topic('all_eatlas_experiment_metadata').collect() )
        .mix( Channel.topic('filtered_eatlas_experiment_metadata').collect() )
        .set { ch_multiqc_files }

    MULTIQC_WORKFLOW( ch_multiqc_files )

    MULTIQC_WORKFLOW.out.report.toList().set { multiqc_report }

    // -----------------------------------------------------------------
    // DASH APPLICATION
    // -----------------------------------------------------------------

    DASH_APP(
        ch_all_counts,
        ch_whole_design,
        ch_all_genes_statistics
    )


    emit:
        multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
