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
include { MERGE_COMPUTE_STATS                    } from '../subworkflows/local/merge_compute_stats'
include { MULTIQC_WORKFLOW                       } from '../subworkflows/local/multiqc'



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
        .view()
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
        // MERGE DATA AND COMPUTE VARIOUS STATISTICS
        // -----------------------------------------------------------------

        MERGE_COMPUTE_STATS (
            DATA_CLEANSING.out.cleaned_counts,
            IDMAPPING.out.gene_metadata,
            IDMAPPING.out.gene_id_mapping

        )

        MERGE_COMPUTE_STATS.out.top_stable_genes_summary.set { ch_top_stable_genes_summary }
        MERGE_COMPUTE_STATS.out.all_genes_statistics.set { ch_all_genes_statistics }
        MERGE_COMPUTE_STATS.out.top_stable_genes_transposed_counts.set { ch_top_stable_genes_transposed_counts }

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


    emit:
        multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
