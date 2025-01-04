//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

import org.yaml.snakeyaml.Yaml

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'
include { workflowVersionToYAML     } from '../../nf-core/utils_nfcore_pipeline'
include { samplesheetToList         } from 'plugin/nf-schema'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs)
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}



/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Setting these functions as temporary replacements of the native softwareVersionsToYAML / processVersionsFromYAML

//
// Get software versions for pipeline
//
def customProcessVersionsFromYAML(yaml_file) {
    Yaml yaml = new Yaml()
    versions = yaml.load(yaml_file)
    return yaml.dumpAsMap(versions).trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def customSoftwareVersionsToYAML(versions) {
    return Channel.of(workflowVersionToYAML())
            .concat(
                versions
                .unique()
                .map {
                    name, tool, version -> [ name.tokenize(':').last(), [ tool, version ] ]
                }
                .groupTuple()
                .map {
                    processName, toolInfo ->
                        def toolVersions = toolInfo.collect { tool, version -> "    ${tool}: ${version}" }.join('\n')
                        "${processName}:\n${toolVersions}\n"
                }
                .map { customProcessVersionsFromYAML(it) }
            )
}

//
// Parses files from input dataset and creates two subchannels raw and normalized
// with elements like [meta, count_file, normalised]
//
def parseInputDatasets(samplesheet) {
    return Channel.fromList( samplesheetToList(samplesheet, "assets/schema_input.json") )
            .map {
                item ->
                    def (count_file, design_file, normalised) = item
                    meta = [dataset: count_file.getBaseName(), design: design_file]
                    [meta, count_file, normalised]
            }
            .branch {
                item ->
                    normalised: item[2] == true
                    raw: item[2] == false
            }
}

//
// Get Expression Atlas Batch ID (accession + data_type) from file stem
//
def augmentWithDatasetId( ch_files ) {
    return ch_files
            .map {
                file ->
                    [[dataset: file.getSimpleName()], file]
            }
}

//
// Groups design and data files by accession and data_type
//
def groupFilesByDatasetId(ch_design, ch_counts) {
    return ch_design
        .concat( ch_counts ) // puts counts at the end of the resulting channel
        .groupTuple() // groups by dataset ID; design files are necessarily BEFORE count files
        .filter {
            it.get(1).size() == 2 // only groups with two files
        }
        .filter { // only groups with first file as design file and second one as count file
            meta, files ->
                files.get(0).name.endsWith('.design.csv') && !files.get(1).name.endsWith('.design.csv')
        }
        .map { // putting design file in meta
            meta, files ->
                def newMeta = meta + [design: files[0]]
                [newMeta, files[1]]
        }
}
