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


include { parseInputDatasets                     } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { customSoftwareVersionsToYAML           } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { validateInputParameters                } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                       } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STABLEEXPRESSION {

    take:
    ch_input_datasets


    main:

    ch_multiqc_files = Channel.empty()

    species = params.species.split(' ').join('_')

    // -----------------------------------------------------------------
    // FETCH EXPRESSION ATLAS DATASETS IF NEEDED
    // -----------------------------------------------------------------

    EXPRESSIONATLAS_FETCHDATA( species )

    // putting all datasets together (local datasets + Expression Atlas datasets)
    ch_input_datasets
        .concat( EXPRESSIONATLAS_FETCHDATA.out.downloaded_datasets )
        .set { ch_datasets }

    // -----------------------------------------------------------------
    // IDMAPPING
    // -----------------------------------------------------------------

    IDMAPPING ( ch_datasets, species )

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
    // MODULE: Merge count files and design files and filter out zero counts
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
    // MODULE: Gene statistics
    // -----------------------------------------------------------------

    GENE_STATISTICS(
        MERGE_DATA.out.all_counts,
        IDMAPPING.out.gene_metadata.collect(),
        IDMAPPING.out.gene_id_mapping.collect(),
        params.nb_top_gene_candidates,
        ch_ks_stats,
        params.ks_pvalue_threshold
    )

    ch_multiqc_files = ch_multiqc_files
                        .mix( GENE_STATISTICS.out.top_stable_genes_summary.collect() )
                        .mix( GENE_STATISTICS.out.all_statistics.collect() )
                        .mix( GENE_STATISTICS.out.top_stable_genes_transposed_counts.collect() )
                        .mix( MERGE_DATA.out.gene_count_statistics.collect() )
                        .mix( MERGE_DATA.out.skewness_statistics.collect() )
                        .mix( ch_ks_stats.collect() )
                        .mix ( MERGE_DATA.out.distribution_correlations.collect() )

    MULTIQC_WORKFLOW( ch_multiqc_files )


    emit:
        multiqc_report = MULTIQC_WORKFLOW.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
