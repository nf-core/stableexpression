include { NORMALISATION_RPK        as RPK               } from '../../../modules/local/normalisation/rpk'
include { NORMALISATION_EDGER_LOG2 as EDGER_LOG2        } from '../../../modules/local/normalisation/edger_log2'

/*
========================================================================================
    SUBWORKFLOW TO COMPUTE GeTMM FROM RNA-SEQ COUNTS
========================================================================================
*/

workflow GETMM_LOG2 {

    take:
    ch_datasets
    ch_gene_length_file

    main:

    //
    // first computing RPK from raw counts
    //

    RPK(
        ch_datasets,
        ch_gene_length_file.collect()
    )

    //
    // feeding these RPK counts to edgeR and compute log2
    //

    EDGER_LOG2( RPK.out.counts )


    emit:
    counts = EDGER_LOG2.out.counts

}
