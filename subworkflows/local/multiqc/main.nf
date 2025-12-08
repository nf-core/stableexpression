include { MULTIQC                                } from '../../../modules/nf-core/multiqc'
include { COLLECT_STATISTICS                     } from '../../../modules/local/collect_statistics'

include { methodsDescriptionText                 } from '../utils_nfcore_stableexpression_pipeline'
include { paramsSummaryMultiqc                   } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                       } from 'plugin/nf-schema'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow MULTIQC_WORKFLOW {

    take:
    ch_multiqc_files
    ch_versions

    main:

    // ------------------------------------------------------------------------------------
    // STATS
    // ------------------------------------------------------------------------------------

    channel.topic('id_mapping_stats')
        .collectFile(
            name: 'id_mapping_stats.csv',
            seed: "Dataset,Nb mapped,Nb unmapped",
            newLine: true,
            storeDir: "${params.outdir}/statistics/"
        ) {
            item -> "${item[0]},${item[1]},${item[2]}"
        }
        .set { ch_id_mapping_stats }

    channel.topic('skewness')
        .map { dataset, file -> "${dataset},${file.readLines()[0]}" } // concatenate dataset name with skewness values
        .collectFile(
            name: 'skewness.csv',
            newLine: true,
            storeDir: "${params.outdir}/statistics/"
        )
        .set { ch_skewness }

    channel.topic('ratio_zeros')
        .map { dataset, file -> "${dataset},${file.readLines()[0]}" } // concatenate dataset name with skewness values
        .collectFile(
            name: 'ratio_zeros.csv',
            newLine: true,
            storeDir: "${params.outdir}/statistics/"
        )
        .set { ch_ratio_zeros }

    ch_to_collect = ch_skewness.mix( ch_ratio_zeros )
    COLLECT_STATISTICS( ch_to_collect )

    // ------------------------------------------------------------------------------------
    // FAILURE / WARNING REPORTS
    // ------------------------------------------------------------------------------------

    channel.topic('eatlas_failure_reason')
        .map { accession, file -> [ accession, file.readLines()[0] ] }
        .collectFile(
            name: 'eatlas_failure_reasons.csv',
            seed: "Accession,Reason",
            newLine: true,
            storeDir: "${params.outdir}/errors/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_eatlas_failure_reasons }

    channel.topic('eatlas_warning_reason')
        .map { accession, file -> [ accession, file.readLines()[0] ] }
        .collectFile(
            name: 'eatlas_warning_reasons.csv',
            seed: "Accession,Reason",
            newLine: true,
            storeDir: "${params.outdir}/warnings/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_eatlas_warning_reasons }

    channel.topic('geo_failure_reason')
        .map { accession, file -> [ accession, file.readLines()[0] ] }
        .collectFile(
            name: 'geo_failure_reasons.csv',
            seed: "Accession,Reason",
            newLine: true,
            storeDir: "${params.outdir}/errors/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_geo_failure_reasons }

    channel.topic('geo_warning_reason')
        .map { accession, file -> [ accession, file.readLines()[0] ] }
        .collectFile(
            name: 'geo_warning_reasons.csv',
            seed: "Accession,Reason",
            newLine: true,
            storeDir: "${params.outdir}/warnings/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_geo_warning_reasons }

    channel.topic('id_cleaning_failure_reason')
        .map { dataset, file -> [ dataset, file.readLines()[0] ] }
        .collectFile(
            name: 'id_cleaning_failure_reasons.tsv',
            seed: "Dataset\tReason",
            newLine: true,
            storeDir: "${params.outdir}/errors/"
        ) {
            item -> "${item[0]}\t${item[1]}"
        }
        .set { ch_id_cleaning_failure_reasons }

    channel.topic('renaming_warning_reason')
        .map { dataset, file -> [ dataset, file.readLines()[0] ] }
        .collectFile(
            name: 'renaming_warning_reasons.tsv',
            seed: "Dataset\tReason",
            newLine: true,
            storeDir: "${params.outdir}/warnings/"
        ) {
            item -> "${item[0]}\t${item[1]}"
        }
        .set { ch_id_mapping_warning_reasons }

    channel.topic('renaming_failure_reason')
        .map { dataset, file -> [ dataset, file.readLines()[0] ] }
        .collectFile(
            name: 'renaming_failure_reasons.tsv',
            seed: "Dataset\tReason",
            newLine: true,
            storeDir: "${params.outdir}/errors/"
        ) {
            item -> "${item[0]}\t${item[1]}"
        }
        .set { ch_id_mapping_failure_reasons }

    channel.topic('normalisation_warning_reason')
        .map { dataset, file -> [ dataset, file.readLines()[0] ] }
        .collectFile(
            name: 'normalisation_warning_reasons.tsv',
            seed: "Dataset\tReason",
            newLine: true,
            storeDir: "${params.outdir}/warnings/"
        ) {
            item -> "${item[0]}\t${item[1]}"
        }
        .set { ch_normalisation_warning_reasons }

    channel.topic('normalisation_failure_reason')
        .map { dataset, file -> [ dataset, file.readLines()[0] ] }
        .collectFile(
            name: 'normalisation_failure_reasons.tsv',
            seed: "Dataset\tReason",
            newLine: true,
            storeDir: "${params.outdir}/errors/"
        ) {
            item -> "${item[0]}\t${item[1]}"
        }
        .set { ch_normalisation_failure_reasons }


    // ------------------------------------------------------------------------------------
    // MULTIQC FILES
    // ------------------------------------------------------------------------------------

    ch_multiqc_files
        .mix( channel.topic('eatlas_all_datasets').collect() )
        .mix( channel.topic('eatlas_selected_datasets').collect() )
        .mix( channel.topic('geo_all_datasets').collect() )
        .mix( channel.topic('geo_selected_datasets').collect() )
        .mix( channel.topic('geo_rejected_datasets').collect() )
        .mix( COLLECT_STATISTICS.out.csv )
        .mix( ch_id_mapping_stats )
        .mix( ch_eatlas_failure_reasons )
        .mix( ch_eatlas_warning_reasons )
        .mix( ch_geo_failure_reasons )
        .mix( ch_geo_warning_reasons )
        .mix( ch_id_cleaning_failure_reasons )
        .mix( ch_id_mapping_warning_reasons )
        .mix( ch_id_mapping_failure_reasons )
        .mix( ch_normalisation_failure_reasons )
        .mix( ch_normalisation_warning_reasons )
        .set { ch_multiqc_files }

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

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'stableexpression_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    // ------------------------------------------------------------------------------------
    // CONFIG
    // ------------------------------------------------------------------------------------

    summary_params = paramsSummaryMap( workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    channel.value( methodsDescriptionText( ch_multiqc_custom_methods_description ) )
        .collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
        .set { ch_methods_description_file }

    ch_multiqc_files
        .mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )
        .mix( ch_collated_versions )
        .mix( ch_methods_description_file )
        .set { ch_multiqc_files }

    ch_multiqc_config = channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    report = MULTIQC.out.report
}
