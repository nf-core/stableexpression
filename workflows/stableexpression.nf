nextflow.enable.dsl = 2
nextflow.preview.topic = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_GETACCESSIONS          } from '../modules/local/expressionatlas/getaccessions/main'
include { EXPRESSIONATLAS_GETDATA                } from '../modules/local/expressionatlas/getdata/main'
include { DESEQ2_NORMALIZE                       } from '../modules/local/deseq2/normalize/main'
include { EDGER_NORMALIZE                        } from '../modules/local/edger/normalize/main'
include { GPROFILER_IDMAPPING                    } from '../modules/local/gprofiler/idmapping/main'
include { VARIATION_COEFFICIENT                  } from '../modules/local/variation_coefficient/main'

include { customSoftwareVersionsToYAML           } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { paramsSummaryMap                       } from 'plugin/nf-schema'
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

    //
    // Initializing channels
    //

    def species = params.species.split(' ').join('_')
    ch_species = Channel.value(species)

    ch_normalized_datasets = Channel.empty()
    ch_raw_datasets = Channel.empty()
    ch_accessions = Channel.empty()

    // if input datasets were provided
    if ( params.datasets ) {

        //
        // Parsing input datasets
        //

        // reads list of input datasets from input file
        // and splits them in normalized and raw sub-channels
        Channel.fromList( samplesheetToList(params.datasets, "${projectDir}/assets/schema_input.json") )
            .map {
                item ->
                    def (count_file, design_file, normalized) = item
                    meta = [accession: count_file.name, design: design_file]
                    [meta, count_file, normalized]
            }
            .branch {
                item ->
                    normalized: item[2] == true
                    raw: item[2] == false
            }
            .set { ch_input_datasets }

        // removes the third element ("normalized" column) and adds to the corresponding channel
        ch_normalized_datasets = ch_normalized_datasets.concat(
            ch_input_datasets.normalized.map{ it -> it.take(2) }
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

    // logging accessions if present
    ch_accessions.collect().map { items -> println "Obtained accessions ${items}"}

    //
    // MODULE: Expression Atlas - Get data
    //

    // Downloading Expression Atlas data for each accession in ch_accessions
    EXPRESSIONATLAS_GETDATA( ch_accessions )

    // separating and arranging EXPRESSIONATLAS_GETDATA output in two separate channels (already normalized or raw data)
    ch_normalized_datasets = ch_normalized_datasets.concat(
        EXPRESSIONATLAS_GETDATA.out.normalized.map {
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
    // MODULE: Normalization of raw count datasets (including RNA-seq datasets)
    //

    if ( params.normalization_method == 'deseq2' ) {
        DESEQ2_NORMALIZE(ch_raw_datasets)
        ch_raw_datasets_normalized = DESEQ2_NORMALIZE.out.csv

    } else { // 'edger'
        EDGER_NORMALIZE(ch_raw_datasets)
        ch_raw_datasets_normalized = EDGER_NORMALIZE.out.csv
    }

    // putting all normalized count datasets together
    ch_normalized_datasets.concat( ch_raw_datasets_normalized ).set{ ch_all_normalized }


    //
    // MODULE: ID Mapping
    //

    // tries to map gene IDs to Ensembl IDs whenever possible
    GPROFILER_IDMAPPING( ch_all_normalized.combine(ch_species) )


    //
    // MODULE: Merge count files & compute variation coefficient for each gene
    //

    VARIATION_COEFFICIENT(
        GPROFILER_IDMAPPING.out.renamed.collect(),
        GPROFILER_IDMAPPING.out.metadata.collect(),
        GPROFILER_IDMAPPING.out.mapping.collect()
    )
    ch_output_from_variation_coefficient = VARIATION_COEFFICIENT.out.csv


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
        ch_output_from_variation_coefficient


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
