include { QUANTILE_NORMALISATION                         } from '../../../modules/local/normalisation/quantile'

include { RNASEQ_NORMALISATION                           } from '../rnaseq_normalisation'
include { MICROARRAY_NORMALISATION                       } from '../microarray_normalisation'

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
    quantile_normalisation
    quantile_norm_target_distrib
    gff_file
    gff_url
    gene_length_file

    main:

    // ------------------------------------------------------------------------------------
    // QUANTILE NORMALISATION
    // ------------------------------------------------------------------------------------

    //
    // set all count datasets together on the same common distribution
    // genes are just ranked among the total set of genes, based on expression
    // and assigned a quantile of rank
    // this method is more scalable but definitely less accurate than the per-platform normalisation
    //

    if ( !quantile_normalisation ) {

        QUANTILE_NORMALISATION (
            ch_datasets,
            quantile_norm_target_distrib
        )
        ch_all_normalised = QUANTILE_NORMALISATION.out.counts

    } else {

        ch_datasets = ch_datasets.branch {
            meta, file ->
                rnaseq: meta.platform == 'rnaseq'
                microarray: meta.platform == 'microarray'
            }

        // ------------------------------------------------------------------------------------
        // NORMALISATION OF RNA-SEQ DATA
        // ------------------------------------------------------------------------------------

        RNASEQ_NORMALISATION( ch_datasets.rnaseq )

        // ------------------------------------------------------------------------------------
        // NORMALISATION OF MICROARRAY DATA
        // ------------------------------------------------------------------------------------




    } else {

        ch_all_normalised = ch_normalised_once

    }


    emit:
    normalised          = ch_all_normalised
    normalised_once     = ch_normalised_once
    gene_length_file    = ch_gene_length_file

}
