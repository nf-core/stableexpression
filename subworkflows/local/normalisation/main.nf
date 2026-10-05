include { QUANTILE_NORMALISATION as QUANTILE     } from '../../../modules/local/normalisation/quantile'
include { MERGE_COUNTS as MERGE_BY_PLATFORM      } from '../../../modules/local/merge_counts'
include { MERGE_COUNTS as REMERGE_BY_PLATFORM    } from '../../../modules/local/merge_counts'
include { IMPUTE_MISSING_VALUES                  } from '../../../modules/local/impute_missing_values'
include { RECOMBAT                               } from '../../../modules/local/recombat'
include { SEPARATE_COUNTS                        } from '../../../modules/local/separate_counts'

include { RNASEQ_NORMALISATION as RNASEQ              } from '../rnaseq_normalisation'

include { mergeDesign                                 } from '../utils_nfcore_stableexpression_pipeline'
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
    missing_value_imputer
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

    RNASEQ(
        ch_datasets.rnaseq,
        species,
        ch_valid_gene_ids,
        skip_gene_length_normalisation,
        gff_file,
        gff_url,
        gene_length_file
    )
    ch_normalised_rnaseq_datasets = RNASEQ.out.normalised
    ch_gene_length_file           = RNASEQ.out.gene_length_file

    // -----------------------------------------------------------------
    // MERGE ALL DATASETS TOGETHER FOR EACH PLATFORM SEPARATELY
    // -----------------------------------------------------------------

    ch_rnaseq = ch_normalised_rnaseq_datasets
                    .map { meta, file -> file }
                    .collect( sort: true )
                    .map { files -> [ [ platform: "rnaseq" ], files ] }

    ch_microarray = ch_datasets.microarray
                        .map { meta, file -> file }
                        .collect( sort: true )
                        .map { files -> [ [ platform: "microarray" ], files ] }

    MERGE_BY_PLATFORM( ch_rnaseq.mix( ch_microarray ) )

    ch_counts_merged_by_platform = MERGE_BY_PLATFORM.out.counts

    // -----------------------------------------------------------------
    // IMPUTE MISSING VALUES
    // -----------------------------------------------------------------

    IMPUTE_MISSING_VALUES(
        ch_counts_merged_by_platform,
        missing_value_imputer
    )

    ch_imputed_datasets = IMPUTE_MISSING_VALUES.out.counts

    // -----------------------------------------------------------------
    // MERGE ALL DESIGNS IN A SINGLE TABLE, PLATFORM PER PLATFORM
    // -----------------------------------------------------------------

    ch_rnaseq_whole_design     = mergeDesign(ch_normalised_rnaseq_datasets, "${outdir}/design/", 'rnaseq.whole_design.csv')
    ch_microarray_whole_design = mergeDesign(ch_datasets.microarray,        "${outdir}/design/", 'microarray.whole_design.csv')

    ch_rnaseq_whole_design     = ch_rnaseq_whole_design.map     { file -> [ 'rnaseq' , file ] }
    ch_microarray_whole_design = ch_microarray_whole_design.map { file -> [ 'microarray', file ] }

    // joining merged datasets with their design
    ch_imputed_datasets_with_design = ch_imputed_datasets
                                        .map { meta, file -> [ meta.platform, meta, file ] }
                                        .join( ch_rnaseq_whole_design.mix( ch_microarray_whole_design ) )
                                        .map { platform, meta, file, design -> [ meta, file, design ] }

    if ( !quantile_normalisation ) {

        // -----------------------------------------------------------------
        // RECOMBAT
        // -----------------------------------------------------------------
        // the default is to correct batch effects, platform per platform
        // first, counts are merged together, platform pre platform
        // then, reCombat is applied on the merged dataframe

        RECOMBAT( ch_imputed_datasets_with_design )

        ch_counts_per_platform = RECOMBAT.out.counts

        // -----------------------------------------------------------------
        // LOGGING NB OF EXCLUDED ORPHAN SAMPLES IF > 0
        // -----------------------------------------------------------------

        RECOMBAT.out.nb_excluded_orphan_samples.map { meta, nb_samples ->
            nb_samples = nb_samples.toInteger()
            if ( nb_samples > 0 ) {
                log.warn("${nb_samples} orphan samples were excluded from platform ${meta.platform}")
            }
        }


    } else {

        // ------------------------------------------------------------------------------------
        // SEPARATE COUNTS AGAIN
        // ------------------------------------------------------------------------------------
        // counts were merged previously in order to impute missing VALUES
        // however, here we do not need to have all datasets together to perform quantile normalisation
        // so we choose to separate datasets again into multiple parquet files, each files corresponding to a batch
        // (separating by batch is only pure convenience)

        SEPARATE_COUNTS( ch_imputed_datasets_with_design )

        ch_separated_datasets = SEPARATE_COUNTS.out.counts
                                    .map { meta, file -> [ [ dataset: file.baseName, platform: meta.platform ], file ] }

        // ------------------------------------------------------------------------------------
        // QUANTILE NORMALISATION
        // ------------------------------------------------------------------------------------
        // force all count datasets to have the same common distribution
        // genes are just ranked among the total set of genes, based on expression
        // and assigned a quantile of rank
        // this method is more scalable but definitely less accurate than the per-platform normalisation
        // first, each dataset is quantile normalised independently, then all normalised datasets are merged together

        QUANTILE(
            ch_separated_datasets,
            quantile_norm_target_distrib
        )
        ch_quantile_normalised_counts = QUANTILE.out.counts

        // -----------------------------------------------------------------
        // MERGE ALL DATASETS TOGETHER FOR EACH PLATFORM SEPARATELY (AGAIN)
        // -----------------------------------------------------------------

        ch_qn_rnaseq_counts = ch_quantile_normalised_counts
                                .filter { meta, file -> meta.platform == "rnaseq" }
                                .map { meta, file -> file }
                                .collect( sort: true )
                                .map { files -> [ [ platform: "rnaseq" ], files ] }

        ch_qn_microarray_datasets = ch_quantile_normalised_counts
                                        .filter { meta, file -> meta.platform == "microarray" }
                                        .map { meta, file -> file }
                                        .collect( sort: true )
                                        .map { files -> [ [ platform: "microarray" ], files ] }

        REMERGE_BY_PLATFORM (
            ch_qn_rnaseq_counts.concat( ch_microarray )
        )

        ch_counts_per_platform = REMERGE_BY_PLATFORM.out.counts

    }

    // -----------------------------------------------------------------
    // ASSOCIATE PLATFORM COUNTS WITH THEIR DESIGN
    // -----------------------------------------------------------------

    ch_merged_counts_with_design = ch_counts_per_platform
                                    .map { meta, file -> [ meta.platform, meta, file ] }
                                    .join( ch_rnaseq_whole_design.mix( ch_microarray_whole_design ) )
                                    .map { platform, meta, file, design -> [ meta, file, design ] }


    emit:
    normalised       = ch_merged_counts_with_design
    imputed          = ch_imputed_datasets
    non_imputed      = ch_counts_merged_by_platform
    gene_length_file = ch_gene_length_file

}
