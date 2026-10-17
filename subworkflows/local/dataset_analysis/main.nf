include { COMPUTE_DATASET_STATISTICS as STATISTICS     } from '../../../modules/local/compute_dataset_statistics'
include { FILTER_SAMPLES                               } from '../../../modules/local/filter_samples'


/*
========================================================================================
    SUBWORKFLOW FOR:
        - COMPUTING VARIOUS STATISTICS AT THE DATASET / SAMPLE LEVEL
        - FILTERING UNVALID SAMLES (TOO MANY NULL VALUES / ZEROS)
========================================================================================
*/

workflow DATASET_ANALYSIS {

    take:
    ch_counts
    ch_valid_gene_ids
    max_zero_ratio
    max_null_ratio
    outdir

    main:

    // -----------------------------------------------------------------
    // COMPUTE VARIOUS STATISTICS AT THE SAMPLE LEVEL
    // -----------------------------------------------------------------

    STATISTICS ( ch_counts )

    // -----------------------------------------------------------------
    // FILTER OUT SAMPLES SHOWING TOO MANY NULL VALUES / TOO MANY ZEROS
    // -----------------------------------------------------------------

    FILTER_SAMPLES (
        ch_counts,
        ch_valid_gene_ids,
        max_null_ratio,
        max_zero_ratio
    )

    // -----------------------------------------------------------------
    // GET NUMBER OF NULLS PER SAMPLE
    // -----------------------------------------------------------------

    ch_ratio_nulls_per_sample_file = FILTER_SAMPLES.out.ratio_nulls_per_sample
                                    .splitCsv( header: true )
                                    .map { item -> "${item["sample"]},${item["ratio"]}" }
                                    .collectFile(
                                        name: 'ratio_nulls_per_sample.csv',
                                        seed: "sample,ratio",
                                        newLine: true,
                                        storeDir: "${outdir}/statistics/",
                                        sort: true
                                    )


    emit:
    counts                         = FILTER_SAMPLES.out.counts
    ratio_nulls_per_sample_file    = ch_ratio_nulls_per_sample_file


}
