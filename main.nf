#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/stableexpression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/stableexpression
    Website: https://nf-co.re/stableexpression
    Slack  : https://nfcore.slack.com/channels/stableexpression
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STABLEEXPRESSION        } from './workflows/stableexpression'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_stableexpression_pipeline'

include { getOutputFolder         } from './subworkflows/local/utils_nfcore_stableexpression_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_STABLEEXPRESSION {

    take:
    input_datasets

    main:

    //
    // WORKFLOW: Run pipeline
    //
    STABLEEXPRESSION( input_datasets )

    emit:
    multiqc_report = STABLEEXPRESSION.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.datasets,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_STABLEEXPRESSION (
        PIPELINE_INITIALISATION.out.input_datasets
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_STABLEEXPRESSION.out.multiqc_report
    )

    publish:
        accessions                            = NFCORE_STABLEEXPRESSION.out.accessions
        input                                 = NFCORE_STABLEEXPRESSION.out.input
        downloaded                            = NFCORE_STABLEEXPRESSION.out.downloaded
        id_filtered_renamed                   = NFCORE_STABLEEXPRESSION.out.id_filtered_renamed
        samples_filtered                      = NFCORE_STABLEEXPRESSION.out.samples_filtered
        first_normalisation                   = NFCORE_STABLEEXPRESSION.out.first_normalisation
        quantile_normalised                   = NFCORE_STABLEEXPRESSION.out.quantile_normalised
        annotation                            = NFCORE_STABLEEXPRESSION.out.annotation
        gene_length_file                      = NFCORE_STABLEEXPRESSION.out.gene_length_file
        merged                                = NFCORE_STABLEEXPRESSION.out.merged
        imputed                               = NFCORE_STABLEEXPRESSION.out.imputed
        all_genes_summary                     = NFCORE_STABLEEXPRESSION.out.all_genes_summary
        multiqc_report                        = NFCORE_STABLEEXPRESSION.out.multiqc_report
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {

    accessions {
        path { source, file ->
            file >> "accessions/${source}/"
        }
    }

    input {
        path { meta, file ->
            file        >> getOutputFolder(meta, "0.provided_by_user")
            meta.design >> getOutputFolder(meta, null)
        }
    }

    downloaded {
        path { meta, file ->
            file        >> getOutputFolder(meta, "0.downloaded")
            meta.design >> getOutputFolder(meta, null)
        }
    }

    id_filtered_renamed {
        path { meta, file ->
            file >> getOutputFolder(meta, "1.id_filtered_renamed")
        }
    }

    samples_filtered {
        path { meta, file ->
            file >> getOutputFolder(meta, "2.samples_filtered")
        }
    }

    first_normalisation {
        path { meta, file ->
            file >> getOutputFolder(meta, "3.${params.normalisation_method}_normalised")
        }
    }

    quantile_normalised {
        path { meta, file ->
            file >> getOutputFolder(meta, "4.quantile_normalised")
        }
    }

    annotation {
        path { file ->
            file >> "annotation/"
        }
    }

    gene_length_file {
        path { meta, file ->
            file >> "annotation/"
        }
    }

    merged {
        path { meta, file ->
            file >> "merged_data/"
        }
    }

    imputed {
        path { meta, file ->
            file >> "merged_data/"
        }
    }

    all_genes_summary {
        path { meta, file ->
            file >> "merged_data/"
        }
    }

    multiqc_report {
        path { file ->
            file >> "reporting/"
        }
    }

}
