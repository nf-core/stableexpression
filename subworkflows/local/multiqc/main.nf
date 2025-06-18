include { MULTIQC                                } from '../../../modules/nf-core/multiqc'

include { customSoftwareVersionsToYAML           } from '../utils_nfcore_stableexpression_pipeline'
include { methodsDescriptionText                 } from '../utils_nfcore_stableexpression_pipeline'
include { paramsSummaryMultiqc                   } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                       } from 'plugin/nf-schema'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow MULTIQC_WORKFLOW {

    take:
    ch_multiqc_files

    main:

    //
    // Collate and save software versions
    //

    ch_collated_versions = customSoftwareVersionsToYAML( Channel.topic('versions') )
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'stableexpression_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

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
