/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_GETDATA                } from '../modules/local/expressionatlas/getdata/main'
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
        //
        // MODULE: Run Expression Atlas - Get data
        //
        ch_species = Channel.value(params.species)
        // in case
        if (params.expression_atlas_keywords == null || params.expression_atlas_keywords == []) {
            params.expression_atlas_keywords = [null]
        }
        ch_keywords = Channel.fromList(params.expression_atlas_keywords)
        ch_eatlas_search = ch_species.combine(ch_keywords)

        ch_eatlas_search | EXPRESSIONATLAS_GETDATA

        EXPRESSIONATLAS_GETDATA.out.csv.view()
    }


    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())



    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
