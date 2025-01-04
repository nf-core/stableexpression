nextflow.enable.dsl = 2
nextflow.preview.topic = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_FETCHDATA              } from '../subworkflows/local/expressionatlas_fetchdata/main'

include { DESEQ2_NORMALISE                       } from '../modules/local/deseq2/normalise/main'
include { EDGER_NORMALISE                        } from '../modules/local/edger/normalise/main'
include { GPROFILER_IDMAPPING                    } from '../modules/local/gprofiler/idmapping/main'
include { GENE_VARIATION                         } from '../modules/local/gene_variation/main'

include { parseInputDatasets                     } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { customSoftwareVersionsToYAML           } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { paramsSummaryLog                       } from 'plugin/nf-schema'
include { validateParameters                     } from 'plugin/nf-schema'


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
    //validateParameters()

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

    ch_species = Channel.value( params.species.split(' ').join('_') )

    ch_normalised_datasets = Channel.empty()
    ch_raw_datasets = Channel.empty()

    // if input datasets were provided
    if ( params.datasets ) {

        //
        // Parsing input datasets
        //

        // reads list of input datasets from input file
        // and splits them in normalised and raw sub-channels
        ch_input_datasets = parseInputDatasets( params.datasets )

        ch_normalised_datasets = ch_normalised_datasets.concat(
            ch_input_datasets.normalised.map{ it -> it.take(2) }
        )
        ch_raw_datasets = ch_raw_datasets.concat(
            ch_input_datasets.raw.map{ it -> it.take(2) }
        )
    }

    //
    // SUBWORKFLOW: fetching Expression Atlas datasets if needed
    //

    EXPRESSIONATLAS_FETCHDATA(
        ch_species,
        params.eatlas_accessions,
        params.eatlas_keywords,
        params.fetch_eatlas_accessions
    )

    // putting all normalized and raw datasets together (local datasets + Expression Atlas datasets)
    ch_normalised_datasets = ch_normalised_datasets.concat( EXPRESSIONATLAS_FETCHDATA.out.normalised )
    ch_raw_datasets = ch_raw_datasets.concat( EXPRESSIONATLAS_FETCHDATA.out.raw )

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
