include { FILTER_OUT_SAMPLES_WITH_TOO_MANY_ZEROS          as TOO_MANY_ZEROS                                   } from '../../../modules/local/filter_out_samples/with_too_many_zeros'
include { FILTER_OUT_SAMPLES_WITH_TOO_MANY_MISSING_VALUES as TOO_MANY_MISSING_VALUES                          } from '../../../modules/local/filter_out_samples/with_too_many_missing_values'


/*
========================================================================================
    SUBWORKFLOW TO FILTER OUT UNVALID SAMPLES AND EMIT STATISTICS ABOUT ZEROS / MISSING VALUES
========================================================================================
*/

workflow SAMPLE_FILTERING {

    take:
    ch_counts
    ch_valid_gene_ids
    max_zero_ratio
    max_null_ratio
    outdir

    main:

    // -----------------------------------------------------------------
    // REMOVE SAMPLES WITH TOO MANY ZEROS
    // -----------------------------------------------------------------

    TOO_MANY_ZEROS (
        ch_counts,
        max_zero_ratio
    )

    // -----------------------------------------------------------------
    // REMOVE SAMPLES WITH TOO MANY MISSING VALUES
    // -----------------------------------------------------------------

    TOO_MANY_MISSING_VALUES(
        TOO_MANY_ZEROS.out.counts,
        ch_valid_gene_ids.collect(),
        max_null_ratio
    )

    // -----------------------------------------------------------------
    // GET NUMBER OF NULLS PER SAMPLE
    // -----------------------------------------------------------------

    ch_ratio_nulls_per_sample_file = TOO_MANY_MISSING_VALUES.out.ratio_nulls_per_sample
                                    .splitCsv( header: true )
                                    .collectFile(
                                        name: 'ratio_nulls_per_sample.csv',
                                        seed: "sample,ratio",
                                        newLine: true,
                                        storeDir: "${outdir}/statistics/",
                                        sort: true
                                    )
                                    {
                                        item -> "${item["sample"]},${item["ratio"]}"
                                    }

    emit:
    counts                         = TOO_MANY_MISSING_VALUES.out.counts
    ratio_nulls_per_sample_file    = ch_ratio_nulls_per_sample_file

}
