/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_GETACCESSIONS          } from '../modules/local/expressionatlas/getaccessions/main'
include { EXPRESSIONATLAS_GETDATA                } from '../modules/local/expressionatlas/getdata/main'
include { IDMAPPING                              } from '../modules/local/idmapping/main'
include { MERGE_COUNT_FILES                      } from '../modules/local/merge_count_files/main'
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

    ch_versions = Channel.empty()

    if (params.species) {

        ch_species = Channel.value(params.species)
        // since we combine both channles, we need to ensure that the keywords channel contains at least one element
        if (params.expression_atlas_keywords == null || params.expression_atlas_keywords == []) {
            params.expression_atlas_keywords = [null]
        }
        ch_keywords = Channel.fromList(params.expression_atlas_keywords)

        //
        // MODULE: Expression Atlas - Get accessions
        //

        ch_species.combine(ch_keywords) | EXPRESSIONATLAS_GETACCESSIONS

        //
        // MODULE: Expression Atlas - Get data
        //

        ch_accessions = EXPRESSIONATLAS_GETACCESSIONS.out.csv
                                .splitCsv()
                                .map{ row -> "${row[0]}"}

        EXPRESSIONATLAS_GETDATA(ch_accessions)

        //
        // MODULE: Id mapping
        //

        EXPRESSIONATLAS_GETDATA.out.csv.combine(ch_species) | IDMAPPING

        //
        // MODULE: Run Merge count files
        //

        MERGE_COUNT_FILES(IDMAPPING.out.csv.toList())

        MERGE_COUNT_FILES.out.csv.view()

    }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
