include { RNASEQ_NORMALISATION           } from '../rnaseq_normalisation'
include { CORRECT_BATCH_EFFECTS          } from '../correct_batch_effects'
include { QUANTILE                       } from '../quantile'

include { mergeDesign                    } from '../../../subworkflow/local/utils_nfcore_stableexpression_pipeline'
/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow NORMALISATION {

    take:
    ch_datasets
    species
    ch_valid_gene_ids
    skip_gene_length_normalisation
    quantile_normalisation
    quantile_norm_target_distrib
    gff_file
    gff_url
    gene_length_file
    outdir

    main:


    ch_datasets = ch_datasets.branch {
        meta, file ->
            rnaseq: meta.platform == 'rnaseq'
            microarray: meta.platform == 'microarray'
        }

    // ------------------------------------------------------------------------------------
    // NORMALISATION OF RNA-SEQ DATA
    // MICROARRAY DATA ARE ALREADY NORMALISED
    // ------------------------------------------------------------------------------------

    RNASEQ_NORMALISATION(
        ch_datasets.rnaseq,
        species,
        ch_valid_gene_ids,
        skip_gene_length_normalisation,
        gff_file,
        gff_url,
        gene_length_file
    )
    ch_normalised_rnaseq_datasets = RNASEQ_NORMALISATION.out.counts


    ch_normalised_datasets = ch_normalised_rnaseq_datasets.mix( ch_datasets.microarray )

    // -----------------------------------------------------------------
    // MERGE ALL DESIGNS IN A SINGLE TABLE, PLATFORM PER PLATFORM
    // -----------------------------------------------------------------

    ch_rnaseq_whole_design     = mergeDesign(ch_normalised_rnaseq_datasets, "${outdir}/merged_data/", 'rnaseq.whole_design.csv')
    ch_microarray_whole_design = mergeDesign(ch_datasets.microarray,        "${outdir}/merged_data/", 'microarray.whole_design.csv')

    ch_rnaseq_whole_design     = ch_rnaseq_whole_design.map     { file -> [ 'rnaseq' , file ] }
    ch_microarray_whole_design = ch_microarray_whole_design.map { file -> [ 'microarray', file ] }


    if ( !quantile_normalisation ) {

        // ------------------------------------------------------------------------------------
        // BATCH EFFECT CORRECTION
        // ------------------------------------------------------------------------------------
        // the default is the correct batch effects, platform per platform
        // first, counts are merged together, platform pre platform
        // then, a specific algorithm is applied on the merged dataframe

        CORRECT_BATCH_EFFECTS(
            ch_normalised_datasets,
            ch_rnaseq_whole_design,
            ch_microarray_whole_design
        )
        ch_counts_per_platform = CORRECT_BATCH_EFFECTS.out.corrected_per_platform

    } else {

        // ------------------------------------------------------------------------------------
        // QUANTILE NORMALISATION
        // ------------------------------------------------------------------------------------
        // set all count datasets together on the same common distribution
        // genes are just ranked among the total set of genes, based on expression
        // and assigned a quantile of rank
        // this method is more scalable but definitely less accurate than the per-platform normalisation
        // first, each dataset is quantile normalised independently, then all normalised datasets are merged together

        QUANTILE (
            ch_normalised_datasets,
            ch_rnaseq_whole_design,
            ch_microarray_whole_design,
            quantile_norm_target_distrib
        )
        ch_counts_per_platform = QUANTILE.out.quantiled_normalised_merged_per_platform

    }


    emit:
    counts_per_platform = ch_counts_per_platform

}
