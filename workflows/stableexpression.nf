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

    multiqc_report = Channel.empty()

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
            params.normalisation_method
        )

        EXPRESSION_NORMALISATION.out.normalised_counts.set { ch_normalised_counts }
        EXPRESSION_NORMALISATION.out.dataset_statistics.set { ch_dataset_statistics }

        // -----------------------------------------------------------------
        // MERGE COUNT FILES AND DESIGN FILES AND FILTER OUT ZERO COUNTS
        // -----------------------------------------------------------------

        MERGE_DATA(
            ch_normalised_counts.map {  meta, file -> [file]        }.collect(),
            ch_normalised_counts.map {  meta, file -> [meta.design] }.collect(),
            ch_dataset_statistics.map { meta, file -> [file]        }.collect(),
            params.nb_top_gene_candidates
        )

        MERGE_DATA.out.candidate_gene_counts.set { ch_candidate_gene_counts }
        MERGE_DATA.out.ks_test_statistics.set { ch_ks_stats }

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

        // -----------------------------------------------------------------
        // MULTIQC
        // -----------------------------------------------------------------

        Channel.empty()
            .mix( GENE_STATISTICS.out.top_stable_genes_summary.collect() )
            .mix( GENE_STATISTICS.out.all_statistics.collect() )
            .mix( GENE_STATISTICS.out.top_stable_genes_transposed_counts.collect() )
            .mix( MERGE_DATA.out.gene_count_statistics.collect() )
            .mix( MERGE_DATA.out.skewness_statistics.collect() )
            .mix( ch_ks_stats.collect() )
            .mix ( MERGE_DATA.out.distribution_correlations.collect() )
            .set { ch_multiqc_files }

        MULTIQC_WORKFLOW( ch_multiqc_files )

        MULTIQC_WORKFLOW.out.report.toList().set { multiqc_report }

    }

    emit:
        multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
