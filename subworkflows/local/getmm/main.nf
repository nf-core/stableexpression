include { NORMALISATION_COMPUTE_RPK   as COMPUTE_RPK     } from '../../../modules/local/normalisation/compute_rpk'
include { NORMALISATION_EDGER         as EDGER           } from '../../../modules/local/normalisation/edger'

/*
========================================================================================
    SUBWORKFLOW TO COMPUTE GeTMM FROM RAW COUNTS
========================================================================================
*/

workflow GETMM {

    take:
    ch_datasets
    ch_gene_length_file

    main:

    //
    // first computing RPK from raw counts
    //

    COMPUTE_RPK(
        ch_datasets,
        ch_gene_length_file.collect()
    )

    //
    // feeding these RPK counts to edgeR
    //

    EDGER( COMPUTE_RPK.out.counts )


    emit:
    counts = EDGER.out.counts

}
