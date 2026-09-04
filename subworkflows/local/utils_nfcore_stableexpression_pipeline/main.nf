//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

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
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

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

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/stableexpression ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/stableexpression/blob/main/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/apptainer/singularity/.../institute> --species <species> --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
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
        ch_input_datasets = channel.empty()
    }

    emit:
    input_datasets = ch_input_datasets

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
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

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//


def check_accession(accession) {
    if ( !( accession.startsWith('E-') || accession.startsWith('GSE') ) ) {
        error('Accession ' + accession + ' is not well formated. All accessions should start with "E-" or "GSE".')
    }
}


def check_accession_string(accessions_str) {
    if ( accessions_str != null && accessions_str != "" ) {
        accessions_str.tokenize(',').each { accession ->
            check_accession(accession)
        }
    }
}

def check_accession_file(accession_file) {
    if ( accession_file != null ) {
        def lines = new File(accession_file).readLines()
        lines.each { accession ->
            check_accession(accession)
        }
    }
}

def validateInputParameters(params) {

    // checking that a species has been provided
    if ( !params.species ) {
        error('You must provide a species name')
    }

    // if accessions are provided or excluded, checking that they are well formated
    check_accession_string( params.accessions )
    check_accession_string( params.excluded_accessions )

    check_accession_file( params.accessions_file )
    check_accession_file( params.excluded_accessions_file )

    if ( params.keywords && params.skip_fetch_eatlas_accessions && !params.fetch_geo_accessions ) {
        log.warn "Ignoring keywords as accessions will not be fetched from Expression Atlas or GEO"
    }

    if ( params.gff && params.gff_url ) {
        log.warn "Both gff and gff_url parameters are provided. Using gff."
    }

}

//
// Parses files from input dataset and creates two subchannels raw and normalized
// with elements like [meta, count_file, normalised]
def parseInputDatasets(samplesheet) {
    return channel.fromList( samplesheetToList(samplesheet, "assets/schema_datasets.json") )
            .map {
                item ->
                    def (meta, count_file) = item
                    def new_meta = meta + [dataset: count_file.getBaseName()]
                    [new_meta, count_file]
            }
}


//
// Validate channels from input samplesheet
//

def getHeaderSeparator(header) {
    def separator = header.contains(',') ? "," : header.contains('\t') ? "\t" : " "
    return separator
}

def validateInputSamplesheet( ch_datasets ) {
    // checking that all microarray datasets (if any) are normalised
    ch_datasets
        .filter {
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

    // checking that all count files are well formated (same number of columns in header and rows)
    ch_datasets
        .map { meta, file ->

            def header = ""
            def separator = ""
            def first_row_components = ""

            if (file.extension == 'gz') {
                def gzipInput = new java.util.zip.GZIPInputStream(file.newInputStream())
                def reader = new java.io.BufferedReader(new java.io.InputStreamReader(gzipInput))

                header = reader.readLine()
                separator = getHeaderSeparator(header)
                first_row_components = reader.readLine().split(separator)

                reader.close()
                gzipInput.close()
            } else {
                header = file.withReader { reader -> reader.readLine() }
                separator = getHeaderSeparator(header)

                first_row_components = file.splitCsv( header: false, skip: 1, limit: 1, sep: separator )[0]
            }

            assert header.split(separator).size() == first_row_components.size() : "Header and first row do not have the same number of columns in file ${file}"
        }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def style = 'style="margin-left: 20px;"'
    def citation_text = [
            "Tools used in the workflow included:",
            "<li ${style}>Expression Atlas (Papatheodorou et al. 2018)</li>",
            params["fetch_geo_accessions"] ? "<li ${style}>NCBI GEO (Edgar et al. 2002)</li>" : "",
            params["skip_id_mapping"] ? "" : "<li ${style}>g:Profiler (Reimand et al. 2007)</li>",
            params["normalisation_method"] == "tpm" && !params["gff"] && !params["gene_length"] ? "<li ${style}>Ensembl (Yates et al. 2026)</li>" : "",
            "<li ${style}>Scikit-learn (Pedregosa et al. 2011)</li>",
            "<li ${style}>NormFinder (Andersen et al. 2004)</li>",
            params["skip_genorm"] ? "" : "<li ${style}>GeNorm (Vandesompele et al. 2002)</li>",
            "<li ${style}>MultiQC (Ewels et al. 2016)</li>"
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Papatheodorou, I., Fonseca, N. A., Keays, M., Tang, Y. A., Barrera, E., Bazant, W., Burke, M., Füllgrabe, A., Fuentes, A. M.-P., George, N., Huerta, L., Koskinen, S., Mohammed, S., Geniza, M., Preece, J., Jaiswal, P., Jarnuczak, A. F., Huber, W., Stegle, O., Vizcaino, J. A., Brazma, A., & Petryszak, R. (2018). Expression Atlas: Gene and protein expression across multiple studies and organisms. Nucleic Acids Research, 46(D1), D246–D251. https://doi.org/10.1093/nar/gkx1158</li>",
            params["fetch_geo_accessions"] ? "<li><Edgar, R., Domrachev, M., & Lash, A. E. (2002). Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Research, 30(1), 207–210. https://doi.org/10.1093/nar/30.1.207</li>" : "",
            params["skip_id_mapping"] ? "" : "<li>Reimand, J., Kull, M., Peterson, H., Hansen, J., & Vilo, J. (2007). g:Profiler—A web-based toolset for functional profiling of gene lists from large-scale experiments. Nucleic Acids Research, 35(suppl_2), W193–W200. https://doi.org/10.1093/nar/gkm226</li>",
            params["normalisation_method"] == "tpm" && !params["gff"] && !params["gene_length"] ? "<li>Yates, A. D., Austine-Orimoloye, O., Azov, A. G., Barba, M., Barnes, I., Barrera-Enriquez, V. P., Becker, A., Bennett, R., Berry, A., Bhai, J., Bhurji, S. K., Branco Lins, P. R., Brooks, L., Budhanuru Ramaraju, S., Campbell, L. I., Carbajo Martinez, M., Carpenter, J., Charkhchi, M., Cortes, L. A., … Finn, R. D. (2026). Ensembl 2026. Nucleic Acids Research, 54(D1), D1053–D1060. https://doi.org/10.1093/nar/gkaf1239</li>" : "",
            "<li>Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Prettenhofer, P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., & Duchesnay, É. (2011). Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research, 12(85), 2825–2830.</li>",
            "<li>Andersen, C. L., Jensen, J. L., & Ørntoft, T. F. (2004). Normalization of Real-Time Quantitative Reverse Transcription-PCR Data: A Model-Based Variance Estimation Approach to Identify Genes Suited for Normalization, Applied to Bladder and Colon Cancer Data Sets. Cancer Research, 64(15), 5245–5250. https://doi.org/10.1158/0008-5472.CAN-04-0496</li>",
            params["skip_genorm"] ? "" : "<li>Vandesompele, J., De Preter, K., Pattyn, F., Poppe, B., Van Roy, N., De Paepe, A., & Speleman, F. (2002). Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome Biology, 3(7), research0034.1. https://doi.org/10.1186/gb-2002-3-7-research0034</li>",
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
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
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
                    def meta = [ dataset: file.getSimpleName() ]
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
def augmentMetadata( ch_files ) {
    return ch_files
            .map {
                meta, file ->
                    def norm_state = getNthPartFromEnd(file.name, 3)
                    def normalised = false
                    if ( norm_state == 'normalised' ) {
                        normalised = true
                    } else if ( norm_state == 'raw' ) {
                        normalised = false
                    } else {
                        error("Invalid normalisation state: ${norm_state}")
                    }

                    def platform = getNthPartFromEnd(file.name, 4)
                    def new_meta = meta + [normalised: normalised, platform: platform]
                    [new_meta, file]
            }
}


/*
========================================================================================
    FUNCTIONS FOR CHECKING NB OF DATASETS
========================================================================================
*/

def checkCounts(ch_counts, fetch_geo_accessions) {

    ch_counts.count().map { n ->
        if( n == 0 ) {
            // display a warning if no datasets are found
            def msg_lst = []
            if ( !fetch_geo_accessions ) {
                msg_lst = [
                    "Could not find any readily usable public dataset...",
                    "This might be due to connection issues on the Expression Atlas FTP server.",
                    "If it is the case, please wait for a couple of minutes and run again.",
                    "Alternatively, datasets for your species of interest might not exist on Expression Atlas.",
                    "In this case, you can try to get additional datasets from NCBI GEO Datasets using the --fetch_geo_accessions flag (this feature is still experimental)."
                ]
            } else {
                msg_lst = [
                    "Could not find any readily usable public dataset...",
                    "This might be due to connection issues on the Expression Atlas FTP server.",
                    "If it is the case, please wait for a couple of minutes and run again.",
                    "You can check directly on NCBI GEO Datasets if there are available datasets for this species that you can prepare yourself:",
                    "https://www.ncbi.nlm.nih.gov/gds",
                    "Once you have prepared your own data, you can relaunch the pipeline and provide your prepared count datasets using the --datasets parameter. ",
                    "For more information, see the online documentation at https://nf-co.re/stableexpression."
                ]
            }
            def msg = msg_lst.join("\n").trim()
            error(msg)
        }
    }
}

/*
========================================================================================
    FUNCTION FOR FORMATING OUTPUT FOLDERS
========================================================================================
*/

def getOutputFolder(meta, subfolder) {
    def normalised_status = meta.normalised ? "normalised" : "raw"
    def subfoldername = subfolder ? "${subfolder}/" : ""
    return "datasets/${meta.platform}/${normalised_status}/${meta.dataset}/${subfoldername}"
}
