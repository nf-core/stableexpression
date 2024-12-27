/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_GETACCESSIONS          } from '../modules/local/expressionatlas/getaccessions/main'
include { EXPRESSIONATLAS_GETDATA                } from '../modules/local/expressionatlas/getdata/main'
include { DESEQ2_NORMALIZE                       } from '../modules/local/deseq2/normalize/main'
include { EDGER_NORMALIZE                        } from '../modules/local/edger/normalize/main'
include { IDMAPPING                              } from '../modules/local/gprofiler/idmapping/main'
include { VARIATION_COEFFICIENT                  } from '../modules/local/variation_coefficient/main'

include { paramsSummaryMap                       } from 'plugin/nf-validation'
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

    if (!params.species) {
        error('You must provide a species name')
    }

    if (
        !params.datasets
        && !params.expression_atlas_accessions
        && !params.fetch_expression_atlas_accessions
        ) {
        error('You must provide at least either --datasets or --fetch_expression_atlas_accessions or --expression_atlas_accessions')
    }

    def species = params.species.split(' ').join('_')
    ch_species = Channel.value(species)

    ch_normalized = Channel.empty()
    ch_raw = Channel.empty()
    ch_accessions = Channel.empty()

    if (params.datasets) {

        log.info "Parsing input data"

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
            .set { ch_input }

        // removes the third element ("normalized" column) and adds to the corresponding channel
        ch_normalized = ch_normalized.concat(
            ch_input.normalized.map{ it -> it.take(2) }
        )
        ch_raw = ch_raw.concat(
            ch_input.raw.map{ it -> it.take(2) }
        )

    }

    if (params.expression_atlas_accessions) {

        ch_accessions = Channel.fromList( params.expression_atlas_accessions.tokenize(',') )

    }

    if (params.fetch_expression_atlas_accessions) {

        //
        // MODULE: Expression Atlas - Get accessions
        //

        // keeping the keywords (separated by spaces) as a single string
        ch_keywords = Channel.value( params.expression_atlas_keywords )

        EXPRESSIONATLAS_GETACCESSIONS(ch_species, ch_keywords)

        // appending to accessions provided by the user
        ch_accessions = ch_accessions.concat( EXPRESSIONATLAS_GETACCESSIONS.out.txt.splitText() )

    }

    //
    // MODULE: Expression Atlas - Get data
    //

    EXPRESSIONATLAS_GETDATA(ch_accessions)

    ch_normalized = ch_normalized.concat(
        EXPRESSIONATLAS_GETDATA.out.normalized.map {
            tuple ->
                def (accession, design_file, count_file) = tuple
                meta = [accession: accession, design: design_file]
                [meta, count_file]
        }
    )

    ch_raw = ch_raw.concat(
        EXPRESSIONATLAS_GETDATA.out.raw.map {
            tuple ->
                def (accession, design_file, count_file) = tuple
                meta = [accession: accession, design: design_file]
                [meta, count_file]
            }
    )


    //
    // MODULE: Normalization of raw count datasets (including RNA-seq datasets)
    //

    if (params.normalization_method == 'deseq2') {
        DESEQ2_NORMALIZE(ch_raw)
        ch_raw_normalized = DESEQ2_NORMALIZE.out.csv

    } else {
        EDGER_NORMALIZE(ch_raw)
        ch_raw_normalized = EDGER_NORMALIZE.out.csv
    }

    // putting all normalized count datasets together
    ch_normalized.concat(ch_raw_normalized).set{ ch_all_normalized }

    //
    // MODULE: Id mapping
    //

    IDMAPPING( ch_all_normalized.combine(ch_species) )

    //
    // MODULE: Merge count files & compute variation coefficient for each gene
    //

    VARIATION_COEFFICIENT( IDMAPPING.out.csv.collect() )
    ch_var_coeff = VARIATION_COEFFICIENT.out.csv



    emit:
    ch_var_coeff
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
