include { MULTIQC                                } from '../../../modules/nf-core/multiqc'

include { formatVersionsToYAML                   } from '../utils_nfcore_stableexpression_pipeline'
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
    // FAILURE / WARNING REPORTS
    // ------------------------------------------------------------------------------------

    Channel.topic('eatlas_failure_reason')
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

    Channel.topic('eatlas_warning_reason')
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

    Channel.topic('geo_failure_reason')
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

    Channel.topic('geo_warning_reason')
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



    // ------------------------------------------------------------------------------------
    // MULTIQC FILES
    // ------------------------------------------------------------------------------------

    ch_multiqc_files
        .mix( Channel.topic('eatlas_all_datasets').collect() )
        .mix( Channel.topic('eatlas_selected_datasets').collect() )
        .mix( ch_eatlas_failure_reasons )
        .mix( ch_eatlas_warning_reasons )
        .mix( Channel.topic('geo_all_datasets').collect() )
        .mix( Channel.topic('geo_selected_datasets').collect() )
        .mix( Channel.topic('geo_rejected_datasets').collect() )
        .mix( ch_geo_failure_reasons )
        .mix( ch_geo_warning_reasons )
        .set { ch_multiqc_files }

    // ------------------------------------------------------------------------------------
    // VERSIONS
    // ------------------------------------------------------------------------------------

    // Collate and save software versions obtained from topic channels
    // TODO: use the nf-core functions when they are adapted to channel topics

    // Collate and save software versions
    formatVersionsToYAML ( Channel.topic('versions') )
        .mix ( softwareVersionsToYAML( ch_versions ) ) // mix with versions obtained from emit outputs
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    // ------------------------------------------------------------------------------------
    // CONFIG
    // ------------------------------------------------------------------------------------

    summary_params = paramsSummaryMap( workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    Channel.value( methodsDescriptionText( ch_multiqc_custom_methods_description ) )
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

    ch_multiqc_config = Channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()

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
