include { NORMALISATION_CPM_LOG2 as CPM_LOG2             } from '../../../modules/local/normalisation/cpm_log2'
include { NORMALISATION_TPM_LOG2 as TPM_LOG2             } from '../../../modules/local/normalisation/tpm_log2'
include { QUANTILE_NORMALISATION                         } from '../../../modules/local/normalisation/quantile'

include { GET_TRANSCRIPT_LENGTHS                         } from '../../../subworkflows/local/get_transcript_lengths'
include { GETMM                                          } from '../../../subworkflows/local/getmm'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow EXPRESSION_NORMALISATION {

    take:
    species
    ch_datasets
    ch_valid_gene_ids
    skip_gene_length_normalisation
    skip_quantile_normalisation
    quantile_norm_target_distrib
    gff_file
    gff_url
    gene_length_file

    main:

    ch_gene_length_file = channel.empty()

    ch_datasets = ch_datasets.branch {
        meta, file ->
            rnaseq: meta.platform == 'rnaseq'
            microarray: meta.platform == 'microarray'
        }

    // ------------------------------------------------------------------------------------
    // NORMALISATION OF RNA-SEQ DATA
    // ------------------------------------------------------------------------------------

    ch_rnaseq_datasets = ch_datasets.rnaseq.branch {
        meta, file ->
            raw: meta.normalised == false
            normalised: meta.normalised == true
        }

    if ( !skip_gene_length_normalisation  ) {

        // ------------------------------------------------------------------------------------
        // GET TABLE ASSOCIATING GENE IDS TO THE LENGTH OF THEIR LARGEST TRANSCRIPT
        // ------------------------------------------------------------------------------------

        if ( gene_length_file ) {

            ch_gene_length_file = channel.fromPath( gene_length_file, checkIfExists: true )

        } else {

            // download genome annotation
            // and compute length of the longest transcript gene per gene
            GET_TRANSCRIPT_LENGTHS(
                species,
                ch_valid_gene_ids,
                gff_file,
                gff_url,
                ch_datasets.rnaseq.collect() // used to trigger this subworkflow only if RNA-seq data are present
            )
            ch_gene_length_file = GET_TRANSCRIPT_LENGTHS.out.csv

        }

        // ------------------------------------------------------------------------------------
        // NORMALISATION
        // ------------------------------------------------------------------------------------

        // normalisation on raw counts (preferred option)
        GETMM(
            ch_rnaseq_datasets.raw,
            ch_gene_length_file
        )

        // if some provided counts are already normalised, use TPM instead
        TPM_LOG2( ch_rnaseq_datasets.normalised )

        ch_normalised_rnaseq_datasets = GETMM.out.counts.mix( TPM_LOG2.out.counts )

    } else {
        // Checking that we have either raw or normalised datasets, but not both
        ch_datasets.rnaseq.toSortedList().map { datasets ->
            def normalisation_status_list = datasets.collect{ meta, file -> meta.normalised ? 'normalised': 'raw' }.unique()
            if ( normalisation_status_list.size() != 1 ) {
                error("Skipping normalisation by gene length is allowed only either raw or normalised datasets are supplied, but not both at once.")
            }
        }

        CPM_LOG2( ch_datasets.rnaseq )
        ch_normalised_rnaseq_datasets = CPM_LOG2.out.counts

    }


    ch_normalised_once = ch_datasets.microarray.mix( ch_normalised_rnaseq_datasets )


    // ------------------------------------------------------------------------------------
    // QUANTILE NORMALISATION
    // ------------------------------------------------------------------------------------

    //
    // set all normalised count datasets together on the common distribution
    //

    if ( !skip_quantile_normalisation ) {

        QUANTILE_NORMALISATION (
            ch_normalised_once,
            quantile_norm_target_distrib
        )
        ch_all_normalised = QUANTILE_NORMALISATION.out.counts

    } else {

        ch_all_normalised = ch_normalised_once

    }


    emit:
    normalised          = ch_all_normalised
    normalised_once     = ch_normalised_once
    gene_length_file    = ch_gene_length_file

}
