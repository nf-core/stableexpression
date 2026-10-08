include { NORMALISATION_TPM_LOG2 as TPM_LOG2             } from '../../../modules/local/normalisation/tpm_log2'
include { QUANTILE_NORMALISATION                         } from '../../../modules/local/normalisation/quantile'
include { NORMALISATION_EDGER_LOG2 as EDGER_LOG2        } from '../../../modules/local/normalisation/edger_log2'

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

    ch_branched_rnaseq_datasets = ch_rnaseq_datasets.branch {
        meta, file ->
            raw: meta.normalised == false
            normalised: meta.normalised == true
        }

    ch_annotation       = channel.empty()
    ch_gene_length_file = channel.empty()

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
            ch_annotation       = GET_TRANSCRIPT_LENGTHS.out.annotation

        }

        // ------------------------------------------------------------------------------------
        // NORMALISATION OR RAW COUNT DATASETS USING GENE LENGTH
        // ------------------------------------------------------------------------------------

        // normalisation on raw counts (preferred option)
        GETMM_LOG2(
            ch_branched_rnaseq_datasets.raw,
            ch_gene_length_file
        )

        ch_normalised_rnaseq_raw_datasets = GETMM_LOG2.out.counts

    } else {

        // ------------------------------------------------------------------------------------
        // NORMALISATION OR RAW COUNT DATASETS WITHOUT GENE LENGTH
        // THIS INTRODUCES A BIAS DUE TO GENE LENGTH, WHICH WILL HAVE AN IMPACT ON THE EXPRESSION LEVEL OF EACH GENE
        // HOWEVER, IT SHOULD NOT HAVE A MAJOR IMPACT ON STABILITY ASSESSMENT
        // ------------------------------------------------------------------------------------

        EDGER_LOG2( ch_branched_rnaseq_datasets.raw )

        ch_normalised_rnaseq_raw_datasets = EDGER_LOG2.out.counts

    }

    // ------------------------------------------------------------------------------------
    // DATASETS ALREADY NORMALISED ARE MAPPED TO TPM (WHEN POSSIBLE)
    // AND LOG2 + 1 IS COMPUTED ON THEM
    // ------------------------------------------------------------------------------------

    // if some provided counts are already normalised, use TPM instead
    TPM_LOG2( ch_branched_rnaseq_datasets.normalised )

    // ------------------------------------------------------------------------------------
    // PUTTING ALL RNASEQ DATASETS TOGETHER
    // ------------------------------------------------------------------------------------

    ch_normalised_rnaseq_datasets = ch_normalised_rnaseq_raw_datasets.mix( TPM_LOG2.out.counts )


    emit:
    normalised          = ch_normalised_rnaseq_datasets
    annotation          = ch_annotation
    gene_length_file    = ch_gene_length_file

}
