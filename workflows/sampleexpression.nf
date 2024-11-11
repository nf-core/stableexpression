/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_GETACCESSIONS          } from '../modules/local/expressionatlas/getaccessions/main'
include { EXPRESSIONATLAS_GETDATA                } from '../modules/local/expressionatlas/getdata/main'
include { DESEQ2_NORMALIZE                       } from '../modules/local/deseq2/normalize/main'
include { EDGER_NORMALIZE                        } from '../modules/local/edger/normalize/main'
include { IDMAPPING                              } from '../modules/local/idmapping/main'
include { MERGE_COUNT_FILES                      } from '../modules/local/merge_count_files/main'
include { MERGE_DESIGNS                          } from '../modules/local/merge_designs/main'
include { VARIATION_COEFFICIENT                  } from '../modules/local/variation_coefficient/main'
include { paramsSummaryMap                       } from 'plugin/nf-validation'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMPLEEXPRESSION {

    //
    // Checking input parameters
    //

    if (!params.species && params.expression_atlas_keywords) {
        error('You must provide a species name if you specify expression atlas keywords')
    }

    if (!params.count_dataset && !params.species) {
        error('You must provide at least either count datasets or a species name')
    }

    if (params.species) {

        def species = params.species.split(' ').join('_')

        ch_species = Channel.value(species)
        ch_keywords = Channel.value(params.expression_atlas_keywords)

        //
        // MODULE: Expression Atlas - Get accessions
        //

        EXPRESSIONATLAS_GETACCESSIONS(ch_species, ch_keywords)

        //
        // MODULE: Expression Atlas - Get data
        //

        ch_accessions = EXPRESSIONATLAS_GETACCESSIONS
                            .out
                            .csv
                            .splitCsv()
                            .map{ row -> "${row[0]}"}

        EXPRESSIONATLAS_GETDATA(ch_accessions)
        EXPRESSIONATLAS_GETDATA.out.normalized.set { ch_normalized }

        //
        // MODULE: Merge the design CSV files into one single file
        //

        MERGE_DESIGNS(EXPRESSIONATLAS_GETDATA.out.design.collect())
        ch_design = MERGE_DESIGNS.out.csv

        //
        // MODULE: Normalization of raw count datasets (including RNA-seq datasets)
        //

        if (params.normalization_method == 'deseq2') {
            DESEQ2_NORMALIZE(
                EXPRESSIONATLAS_GETDATA.out.raw,
                ch_design
            )
            ch_raw_normalized = DESEQ2_NORMALIZE.out.csv

        } else {
            EDGER_NORMALIZE(
                EXPRESSIONATLAS_GETDATA.out.raw,
                ch_design
            )
            ch_raw_normalized = EDGER_NORMALIZE.out.csv
        }

        // putting all normalized count datasets together
        ch_all_normalized = ch_normalized.concat(ch_raw_normalized).map { it -> it[1]}
        ch_all_normalized.view()

        //
        // MODULE: Id mapping
        //

        ch_all_normalized.combine(ch_species) | IDMAPPING

        //
        // MODULE: Run Merge count files
        //

        MERGE_COUNT_FILES(IDMAPPING.out.csv.collect())

        //
        // MODULE: Compute variation coefficient for each gene
        //

        VARIATION_COEFFICIENT(MERGE_COUNT_FILES.out.csv)
        ch_var_coeff = VARIATION_COEFFICIENT.out.csv

    }




    emit:
    ch_var_coeff
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
