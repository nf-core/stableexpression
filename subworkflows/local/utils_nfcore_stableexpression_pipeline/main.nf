//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

import org.yaml.snakeyaml.Yaml

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { workflowVersionToYAML     } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        params.outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters( params )

    //
    // Create channel from datasets file provided through params.datasets
    //
    if (params.datasets) {
        ch_input_datasets = parseInputDatasets( params.datasets )
        validateInputSamplesheet( ch_input_datasets )
    } else {
        ch_input_datasets = Channel.empty()
    }

    emit:
    input_datasets = ch_input_datasets

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
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
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
//
// Check and validate pipeline parameters
//

def validateInputParameters(params) {

    // checking that a species has been provided
    if ( !params.species ) {
        error('You must provide a species name')
    }

    // if expression atlas accessions are provided, checking that they are well formated
    if ( params.eatlas_accessions ) {
        for ( accession in params.eatlas_accessions.tokenize(',') ) {
            if ( !accession.startsWith('E-') ) {
                error('Expression Atlas accession ' + accession + ' is not well formated. All accessions should start with "E-".')
            }
        }
    }

}

//
// Parses files from input dataset and creates two subchannels raw and normalized
// with elements like [meta, count_file, normalised]
def parseInputDatasets(samplesheet) {
    return Channel.fromList( samplesheetToList(samplesheet, "assets/schema_datasets.json") )
            .map {
                item ->
                    def (meta, count_file) = item
                    new_meta = meta + [dataset: count_file.getBaseName()]
                    [new_meta, count_file]
            }
}


//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    // checking that all microarray datasets (if any) are normalised
    input.filter {
        meta, file ->
            meta.platform == 'microarray' && !meta.normalised
    }
    .count()
    .map { count ->
        if (count > 0) {
            def error_text = [
                "Error: You provided at least one microarray dataset that is not normalised. ",
                "Microarray datasets must already be normalised before being submitted. ",
                "Please perform normalisation (typically using RMA for one-colour intensities / LOESS (limma) for two-colour intensities) and run again."
            ].join(' ').trim()
            error(error_text)
        }
    }
}

//
// Get channel of software versions used in pipeline in YAML format
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Get software versions for pipeline
// temporary replacements of the native processVersionsFromYAML
//
def customProcessVersionsFromYAML(yaml_file) {
    Yaml yaml = new Yaml()
    versions = yaml.load(yaml_file)
    return yaml.dumpAsMap(versions).trim()
}

//
// Get channel of software versions used in pipeline in YAML format
// temporary replacements of the native softwareVersionsToYAML
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



/*
========================================================================================
    FUNCTIONS FOR FORMATTING DATA FETCHED FROM EXPRESSION ATLAS / GEO
========================================================================================
*/

//
// Get Expression Atlas Batch ID (accession + data_type) from file stem
//
def addDatasetIdToMetadata( ch_files ) {
    return ch_files
            .map {
                file ->
                    def meta = [dataset: file.getSimpleName()]
                    [meta, file]
            }
}

//
// Groups design and data files by accession and data_type
// Design and count files have necessarily the same dataset ID (same file stem)
//
def groupFilesByDatasetId(ch_design, ch_counts) {
    return ch_design
        .concat( ch_counts ) // puts counts at the end of the resulting channel
        .groupTuple() // groups by dataset ID; design files are necessarily BEFORE count files
        .filter {
            it.get(1).size() == 2 // only groups with two files
        }
        .filter { // only groups with first file as design file and second one as count fileWARN: java.net.ConnectException: Connexion refusée
            meta, files ->
                files.get(0).name.endsWith('.design.csv') && !files.get(1).name.endsWith('.design.csv')
        }
        .map { // putting design file in meta
            meta, files ->
                def new_meta = meta + [design: files[0]]
                [new_meta, files[1]]
        }
}

def getNthPartFromEnd(String s, int n) {
    def tokens = s.tokenize('.')
    return tokens[tokens.size() - n]
}

//
// Add normalised: true / false in meta
//
def augmentToMetadata( ch_files ) {
    return ch_files
            .map {
                meta, file ->
                    if ( getNthPartFromEnd(file.name, 3) == 'raw' ) {
                        meta.normalised = false
                    } else {
                        meta.normalised = true
                    }
                    meta.platform = getNthPartFromEnd(file.name, 4)
                    [meta, file]
            }
}




