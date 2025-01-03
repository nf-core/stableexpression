nextflow.enable.dsl = 2
nextflow.preview.topic = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_GETACCESSIONS          } from '../modules/local/expressionatlas/getaccessions/main'
include { EXPRESSIONATLAS_GETDATA                } from '../modules/local/expressionatlas/getdata/main'
include { DESEQ2_NORMALISE                       } from '../modules/local/deseq2/normalise/main'
include { EDGER_NORMALISE                        } from '../modules/local/edger/normalise/main'
include { GPROFILER_IDMAPPING                    } from '../modules/local/gprofiler/idmapping/main'
include { GENE_VARIATION                         } from '../modules/local/gene_variation/main'

include { customSoftwareVersionsToYAML           } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { paramsSummaryLog                       } from 'plugin/nf-schema'
include { validateParameters                     } from 'plugin/nf-schema'
include { samplesheetToList                      } from 'plugin/nf-schema'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STABLEEXPRESSION {

    //
    // Checking input parameters
    //

    // Validate input parameters
    validateParameters()

    if ( !params.species ) {
        error('You must provide a species name')
    }

    if (
        !params.datasets
        && !params.eatlas_accessions
        && !params.fetch_eatlas_accessions
        && !params.eatlas_keywords
        ) {
        error('You must provide at least either --datasets or --fetch_eatlas_accessions or --eatlas_accessions or --eatlas_keywords')
    }

    // Print summary of supplied parameters
    // log.info paramsSummaryLog(workflow)

    //
    // Initializing channels
    //

    def species = params.species.split(' ').join('_')
    ch_species = Channel.value( species )

    ch_normalised_datasets = Channel.empty()
    ch_raw_datasets = Channel.empty()
    ch_accessions = Channel.empty()

    // if input datasets were provided
    if ( params.datasets ) {

        //
        // Parsing input datasets
        //

        // reads list of input datasets from input file
        // and splits them in normalised and raw sub-channels
        Channel.fromList(samplesheetToList(params.datasets, "assets/schema_input.json"))
            .map {
                item ->
                    def (count_file, design_file, normalised) = item
                    meta = [accession: count_file.name, design: design_file]
                    [meta, count_file, normalised]
            }
            .branch {
                item ->
                    normalised: item[2] == true
                    raw: item[2] == false
            }
            .set { ch_input_datasets }

        // removes the third element ("normalised" column) and adds to the corresponding channel
        ch_normalised_datasets = ch_normalised_datasets.concat(
            ch_input_datasets.normalised.map{ it -> it.take(2) }
        )
        ch_raw_datasets = ch_raw_datasets.concat(
            ch_input_datasets.raw.map{ it -> it.take(2) }
        )

    }

    // parsing Expression Atlas accessions if provided
    if ( params.eatlas_accessions ) {

        // parsing accessions from provided parameter
        ch_accessions = Channel.fromList( params.eatlas_accessions.tokenize(',') )

    }


    // fetching Expression Atlas accessions if applicable
    if ( params.fetch_eatlas_accessions || params.eatlas_keywords ) {

        //
        // MODULE: Expression Atlas - Get accessions
        //

        // keeping the keywords (separated by spaces) as a single string
        ch_keywords = Channel.value( params.eatlas_keywords )

        // getting Expression Atlas accessions given a species name and keywords
        // keywords can be an empty string
        EXPRESSIONATLAS_GETACCESSIONS( ch_species, ch_keywords )

        // appending to accessions provided by the user
        // ensures that no accessions is present twice (provided by the user and fetched from E. Atlas)
        ch_accessions = ch_accessions
                            .concat( EXPRESSIONATLAS_GETACCESSIONS.out.txt.splitText() )
                            .unique()

    }

    //
    // MODULE: Expression Atlas - Get data
    //

    // Downloading Expression Atlas data for each accession in ch_accessions
    EXPRESSIONATLAS_GETDATA( ch_accessions )

    // separating and arranging EXPRESSIONATLAS_GETDATA output in two separate channels (already normalised or raw data)
    ch_normalised_datasets = ch_normalised_datasets.concat(
        EXPRESSIONATLAS_GETDATA.out.normalised.map {
            accession, design_file, count_file ->
                meta = [accession: accession, design: design_file]
                [meta, count_file]
        }
    )

    ch_raw_datasets = ch_raw_datasets.concat(
        EXPRESSIONATLAS_GETDATA.out.raw.map {
            accession, design_file, count_file ->
                meta = [accession: accession, design: design_file]
                [meta, count_file]
            }
    )


    //
    // MODULE: normalisation of raw count datasets (including RNA-seq datasets)
    //

    if ( params.normalisation_method == 'deseq2' ) {
        DESEQ2_NORMALISE( ch_raw_datasets )
        ch_raw_datasets_normalised = DESEQ2_NORMALISE.out.csv

    } else { // 'edger'
        EDGER_NORMALISE( ch_raw_datasets )
        ch_raw_datasets_normalised = EDGER_NORMALISE.out.csv
    }

    // putting all normalised count datasets together
    ch_normalised_datasets.concat( ch_raw_datasets_normalised ).set{ ch_all_normalised }


    //
    // MODULE: ID Mapping
    //

    // tries to map gene IDs to Ensembl IDs whenever possible
    GPROFILER_IDMAPPING( ch_all_normalised.combine(ch_species) )

    //
    // MODULE: Merge count files & compute variation coefficient for each gene
    //

    GENE_VARIATION(
        GPROFILER_IDMAPPING.out.renamed.collect(),
        GPROFILER_IDMAPPING.out.metadata.collect(),
        GPROFILER_IDMAPPING.out.mapping.collect(),
        Channel.value( params.gene_variation_method )
    )
    ch_output_from_gene_variation = GENE_VARIATION.out.csv


    //
    // Collate and save software versions
    // TODO: use the nf-core functions when they are adapted to channel topics
    //

    customSoftwareVersionsToYAML( Channel.topic('versions') )
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true
        )

    // only used for nf-test
    emit:
        ch_output_from_gene_variation


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
