include { NORMALISATION_CPM_LOG2 as CPM_LOG2             } from '../../../modules/local/normalisation/cpm_log2'
include { NORMALISATION_TPM_LOG2 as TPM_LOG2             } from '../../../modules/local/normalisation/tpm_log2'
include { QUANTILE_NORMALISATION                         } from '../../../modules/local/normalisation/quantile'

include { GET_TRANSCRIPT_LENGTHS                         } from '../../../subworkflows/local/get_transcript_lengths'
include { GETMM_LOG2                                     } from '../../../subworkflows/local/getmm_log2'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE RNA-SEQ EXPRESSION DATASETS
========================================================================================
*/

workflow RNASEQ_NORMALISATION {

    take:
    ch_rnaseq_datasets
    species
    ch_valid_gene_ids
    skip_gene_length_normalisation
    gff_file
    gff_url
    gene_length_file

    main:

    ch_gene_length_file = channel.empty()

    // ------------------------------------------------------------------------------------
    // AT THIS POINT, ONLY ONE BRANCH (RAW OR NORMALISED) HAS ELEMENTS
    // ------------------------------------------------------------------------------------

    ch_branched_rnaseq_datasets = ch_rnaseq_datasets.branch {
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
            // and compute length of the longest transcript gene per geneRNA-SEQ
            GET_TRANSCRIPT_LENGTHS(
                species,
                ch_valid_gene_ids,
                gff_file,
                gff_url,
                ch_rnaseq_datasets.collect() // used to trigger this subworkflow only if RNA-seq data are present
            )
            ch_gene_length_file = GET_TRANSCRIPT_LENGTHS.out.csv

        }

        // ------------------------------------------------------------------------------------
        // NORMALISATION OR RAW COUNT DATASETS
        // ------------------------------------------------------------------------------------

        // normalisation on raw counts (preferred option)
        GETMM_LOG2(
            ch_branched_rnaseq_datasets.raw,
            ch_gene_length_file
        )

        // ------------------------------------------------------------------------------------
        // RE-NORMALISATION OR ALREADY NORMALISED COUNT DATASETS
        // ------------------------------------------------------------------------------------

        // if some provided counts are already normalised, use TPM instead
        TPM_LOG2( ch_branched_rnaseq_datasets.normalised )

        ch_normalised_rnaseq_datasets = GETMM_LOG2.out.counts.mix( TPM_LOG2.out.counts )

    } else {

        CPM_LOG2( ch_rnaseq_datasets )
        ch_normalised_rnaseq_datasets = CPM_LOG2.out.counts

    }



    emit:
    normalised          = ch_normalised_rnaseq_datasets
    gene_length_file    = ch_gene_length_file

}
