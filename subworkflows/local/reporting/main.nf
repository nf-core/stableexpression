include { AGGREGATE_RESULTS                      } from '../../../modules/local/aggregate_results'
include { DASH_APP                               } from '../../../modules/local/dash_app'
include { COLLECT_STATISTICS                     } from '../../../modules/local/collect_statistics'
include { MULTIQC                                } from '../../../modules/nf-core/multiqc'


include { methodsDescriptionText                 } from '../utils_nfcore_stableexpression_pipeline'
include { paramsSummaryMultiqc                   } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                       } from 'plugin/nf-schema'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow REPORTING {

    take:
    ch_all_counts
    ch_whole_design
    ch_stats_all_genes_with_scores
    ch_platform_statistics
    ch_whole_gene_metadata
    ch_whole_gene_id_mapping
    target_genes
    target_gene_file
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    ch_versions = channel.empty()

    // -----------------------------------------------------------------
    // AGGREGATE ALL RESULTS FOR MULTIQC
    // -----------------------------------------------------------------

    ch_target_gene_file = target_gene_file ? channel.fromPath( target_gene_file, checkIfExists: true ) : channel.empty()

    ch_target_gene_list = channel.fromList( target_genes.tokenize(',') )
                        .mix( ch_target_gene_file.splitText() )
                        .map { it.trim() }
                        .filter { it != "" }
                        .unique()
                        .toSortedList()

    ch_custom_content_multiqc_config_template = channel.fromPath(
                                                    "${projectDir}/assets/multiqc_config.custom_content.template.yaml",
                                                    checkIfExists: true
                                                )

    AGGREGATE_RESULTS (
        ch_all_counts.map{ meta, file -> file }.collect(),
        ch_stats_all_genes_with_scores.collect(),
        ch_platform_statistics.collect(),
        ch_target_gene_list,
        ch_whole_gene_metadata.collect().ifEmpty([]), // handle case where there are no mappings
        ch_whole_gene_id_mapping.collect().ifEmpty([]), // handle case where there are no mappings
        ch_custom_content_multiqc_config_template.collect()
    )

    ch_all_genes_summary                   = AGGREGATE_RESULTS.out.all_genes_summary
    ch_most_stable_genes_summary           = AGGREGATE_RESULTS.out.most_stable_genes_summary
    ch_most_stable_genes_transposed_counts = AGGREGATE_RESULTS.out.most_stable_genes_transposed_counts_filtered
    ch_custom_content_multiqc_config       = AGGREGATE_RESULTS.out.custom_content_multiqc_config

    // -----------------------------------------------------------------
    // DASH APPLICATION
    // -----------------------------------------------------------------

    DASH_APP(
        ch_all_counts.map{ meta, file -> file }.collect(),
        ch_whole_design.collect(),
        ch_all_genes_summary.collect()
    )
    ch_versions = ch_versions.mix ( DASH_APP.out.versions )


    // ------------------------------------------------------------------------------------
    // PREPARING BAR PLOTS
    // ------------------------------------------------------------------------------------

    ch_id_mapping_stats = channel.topic('mqc_id_mapping_stats')
                            .collectFile(
                                name: 'id_mapping_stats.csv',
                                seed: "dataset,final,merged,not_valid,unmapped",
                                newLine: true,
                                storeDir: "${outdir}/statistics/"
                            ) {
                                item -> "${item[0]},${item[1]},${item[2]},${item[3]},${item[4]}"
                            }

    ch_missing_values_filter_stats = channel.topic('mqc_missing_values_filter_stats')
                                        .collectFile(
                                            name: 'missing_values_filter_stats.csv',
                                            seed: "dataset,kept,rejected",
                                            newLine: true,
                                            storeDir: "${outdir}/statistics/"
                                        ) {
                                            item -> "${item[0]},${item[1]},${item[2]}"
                                        }

    ch_zero_values_filter_stats = channel.topic('mqc_zero_values_filter_stats')
                                .collectFile(
                                    name: 'zero_values_filter_stats.csv',
                                    seed: "dataset,kept,rejected",
                                    newLine: true,
                                    storeDir: "${outdir}/statistics/"
                                ) {
                                    item -> "${item[0]},${item[1]},${item[2]}"
                                }

    // ------------------------------------------------------------------------------------
    // PREPARING BOX PLOTS
    // ------------------------------------------------------------------------------------

    ch_skewness         = channel.topic('skewness')
                            .map { dataset, file -> "${dataset},${file.readLines()[0]}" } // concatenate dataset name with skewness values
                            .collectFile(
                                name: 'skewness.csv',
                                newLine: true,
                                sort: true,
                                storeDir: "${outdir}/statistics/"
                            )


    ch_ratio_zeros      = channel.topic('ratio_zeros')
                            .map { dataset, file -> "${dataset},${file.readLines()[0]}" } // concatenate dataset name with ratio values
                            .collectFile(
                                name: 'ratio_zeros.csv',
                                newLine: true,
                                sort: true,
                                storeDir: "${outdir}/statistics/"
                                )

    ch_ratio_nulls      = channel.topic('ratio_nulls')
                            .map { dataset, file -> "${dataset},${file.readLines()[0]}" } // concatenate dataset name with ratio values
                            .collectFile(
                                name: 'ratio_nulls.csv',
                                newLine: true,
                                sort: true,
                                storeDir: "${outdir}/statistics/"
                                )

    ch_stat_files = ch_skewness
                        .mix( ch_ratio_nulls )
                        .mix( ch_ratio_zeros )

    COLLECT_STATISTICS( ch_stat_files )

    // ------------------------------------------------------------------------------------
    // FAILURE / WARNING REPORTS
    // ------------------------------------------------------------------------------------

    ch_eatlas_failure_reasons = channel.topic('eatlas_failure_reason')
                                    .map { accession, file -> [ accession, file.readLines()[0] ] }
                                    .collectFile(
                                        name: 'eatlas_failure_reasons.csv',
                                        seed: "Accession,Reason",
                                        newLine: true,
                                        sort: true,
                                        storeDir: "${outdir}/errors/",
                                    ) {
                                        item -> "${item[0]},${item[1]}"
                                    }

    ch_eatlas_warning_reasons = channel.topic('eatlas_warning_reason')
                                    .map { accession, file -> [ accession, file.readLines()[0] ] }
                                    .collectFile(
                                        name: 'eatlas_warning_reasons.csv',
                                        seed: "Accession,Reason",
                                        newLine: true,
                                        sort: true,
                                        storeDir: "${outdir}/warnings/"
                                    ) {
                                        item -> "${item[0]},${item[1]}"
                                    }

    ch_geo_failure_reasons = channel.topic('geo_failure_reason')
                                .map { accession, file -> [ accession, file.readLines()[0] ] }
                                .collectFile(
                                    name: 'geo_failure_reasons.csv',
                                    seed: "Accession,Reason",
                                    newLine: true,
                                    sort: true,
                                    storeDir: "${outdir}/errors/"
                                ) {
                                    item -> "${item[0]},${item[1]}"
                                }


    ch_geo_warning_reasons = channel.topic('geo_warning_reason')
                                .map { accession, file -> [ accession, file.readLines()[0] ] }
                                .collectFile(
                                    name: 'geo_warning_reasons.csv',
                                    seed: "Accession,Reason",
                                    newLine: true,
                                    sort: true,
                                    storeDir: "${outdir}/warnings/"
                                ) {
                                    item -> "${item[0]},${item[1]}"
                                }

    ch_id_cleaning_failure_reasons = channel.topic('id_cleaning_failure_reason')
                                        .map { dataset, file -> [ dataset, file.readLines()[0] ] }
                                        .collectFile(
                                            name: 'id_cleaning_failure_reasons.tsv',
                                            seed: "Dataset\tReason",
                                            newLine: true,
                                            sort: true,
                                            storeDir: "${outdir}/errors/"
                                        ) {
                                            item -> "${item[0]}\t${item[1]}"
                                        }

    ch_id_mapping_warning_reasons = channel.topic('renaming_warning_reason')
                                        .map { dataset, file -> [ dataset, file.readLines()[0] ] }
                                        .collectFile(
                                            name: 'renaming_warning_reasons.tsv',
                                            seed: "Dataset\tReason",
                                            newLine: true,
                                            sort: true,
                                            storeDir: "${outdir}/warnings/"
                                        ) {
                                            item -> "${item[0]}\t${item[1]}"
                                        }

    ch_id_mapping_failure_reasons = channel.topic('renaming_failure_reason')
                                        .map { dataset, file -> [ dataset, file.readLines()[0] ] }
                                        .collectFile(
                                            name: 'renaming_failure_reasons.tsv',
                                            seed: "Dataset\tReason",
                                            newLine: true,
                                            sort: true,
                                            storeDir: "${outdir}/errors/"
                                        ) {
                                            item -> "${item[0]}\t${item[1]}"
                                        }

    ch_normalisation_warning_reasons = channel.topic('normalisation_warning_reason')
                                            .map { dataset, file -> [ dataset, file.readLines()[0] ] }
                                            .collectFile(
                                                name: 'normalisation_warning_reasons.tsv',
                                                seed: "Dataset\tReason",
                                                newLine: true,
                                                sort: true,
                                                storeDir: "${outdir}/warnings/"
                                            ) {
                                                item -> "${item[0]}\t${item[1]}"
                                            }

    ch_normalisation_failure_reasons = channel.topic('normalisation_failure_reason')
                                            .map { dataset, file -> [ dataset, file.readLines()[0] ] }
                                            .collectFile(
                                                name: 'normalisation_failure_reasons.tsv',
                                                seed: "Dataset\tReason",
                                                newLine: true,
                                                sort: true,
                                                storeDir: "${outdir}/errors/"
                                            ) {
                                                item -> "${item[0]}\t${item[1]}"
                                            }


    // ------------------------------------------------------------------------------------
    // MULTIQC FILES
    // ------------------------------------------------------------------------------------

    ch_multiqc_files = channel.empty()
                        .mix( ch_most_stable_genes_summary.collect() )                          // single item
                        .mix( ch_all_genes_summary.collect() )                                  // single item
                        .mix( ch_most_stable_genes_transposed_counts.collect() )                // single item
                        .mix( channel.topic('eatlas_all_datasets').toSortedList() )
                        .mix( channel.topic('eatlas_selected_datasets').toSortedList() )
                        .mix( channel.topic('geo_all_datasets').toSortedList() )
                        .mix( channel.topic('geo_selected_datasets').toSortedList() )
                        .mix( channel.topic('geo_rejected_datasets').toSortedList() )
                        .mix( channel.topic('total_gene_id_occurrence_quantiles').toSortedList() )
                        .mix( COLLECT_STATISTICS.out.csv )
                        .mix( ch_id_mapping_stats )
                        .mix( ch_missing_values_filter_stats )
                        .mix( ch_zero_values_filter_stats )
                        .mix( ch_eatlas_failure_reasons )
                        .mix( ch_eatlas_warning_reasons )
                        .mix( ch_geo_failure_reasons )
                        .mix( ch_geo_warning_reasons )
                        .mix( ch_id_cleaning_failure_reasons )
                        .mix( ch_id_mapping_warning_reasons )
                        .mix( ch_id_mapping_failure_reasons )
                        .mix( ch_normalisation_failure_reasons )
                        .mix( ch_normalisation_warning_reasons )


    // ------------------------------------------------------------------------------------
    // VERSIONS
    // ------------------------------------------------------------------------------------

    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
                            .mix(topic_versions_string)
                            .collectFile(
                                storeDir: "${outdir}/pipeline_info",
                                name: 'nf_core_'  +  'stableexpression_software_'  + 'mqc_'  + 'versions.yml',
                                sort: true,
                                newLine: true
                            )

    // ------------------------------------------------------------------------------------
    // PREPARE MULTIQC INPUT
    // ------------------------------------------------------------------------------------

    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    ch_multiqc_custom_config = multiqc_config ?
        channel.fromPath(multiqc_config, checkIfExists: true) :
        channel.empty()

    ch_multiqc_logo          = multiqc_logo ?
        channel.fromPath(multiqc_logo, checkIfExists: true) :
        channel.of([])

    summary_params      = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_files = ch_multiqc_files
        .mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )

    ch_multiqc_custom_methods_description = multiqc_methods_description ?
        file(multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_methods_description     = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    // ------------------------------------------------------------------------------------
    // ADDING KEY TO JOIN ON
    // ------------------------------------------------------------------------------------

    ch_multiqc_file_list = ch_multiqc_files
                            .mix( ch_collated_versions )
                            .mix(
                                ch_methods_description.collectFile(
                                    name: 'methods_description_mqc.yaml',
                                    sort: true
                                )
                            )
                            .flatten()
                            .toSortedList()
                            .map{ list -> [ [id: 'Final report'], list ] }

    ch_multiqc_config_list = ch_multiqc_config
                                .mix( ch_multiqc_custom_config )
                                .mix( ch_custom_content_multiqc_config )
                                .toSortedList()
                                .map{ list -> [ [id: 'Final report'], list ] }

    ch_multiqc_logo = ch_multiqc_logo.map{ file -> [ [id: 'Final report'], file ] }

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_multiqc_input = ch_multiqc_file_list
                        .join( ch_multiqc_config_list )
                        .join( ch_multiqc_logo )
                        .map { meta, files, configs, logo -> [ meta, files, configs, logo , [], [] ] }

    MULTIQC ( ch_multiqc_input )

    emit:
    multiqc_report          = MULTIQC.out.report
    all_genes_summary       = ch_all_genes_summary
}
