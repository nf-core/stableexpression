/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_FETCHDATA              } from '../subworkflows/local/expressionatlas_fetchdata/main'
include { IDMAPPING                              } from '../subworkflows/local/idmapping/main.nf'
include { EXPRESSION_NORMALISATION               } from '../subworkflows/local/expression_normalisation/main.nf'


include { MERGE_DATA                             } from '../modules/local/merge_data/main'
include { GENE_STATISTICS                        } from '../modules/local/gene_statistics/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'

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


    //
    // Collate and save software versions
    // TODO: use the nf-core functions when they are adapted to channel topics
    //

    ch_collated_versions = customSoftwareVersionsToYAML( Channel.topic('versions') )
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'stableexpression_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
        multiqc_report = MULTIQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
