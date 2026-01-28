include { COMPUTE_DATASET_STATISTICS as STATISTICS                     } from '../../../modules/local/compute_dataset_statistics'
include { GET_NB_NULLS_PER_SAMPLE    as NB_MISSING_VALUES              } from '../../../modules/local/get_nb_nulls_per_sample'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow DATASET_ANALYSIS {

    take:
    ch_counts
    ch_valid_gene_ids
    outdir

    main:

    // -----------------------------------------------------------------
    // COMPUTE VARIOUS STATISTICS AT THE SAMPLE LEVEL
    // -----------------------------------------------------------------

    STATISTICS ( ch_counts )

    // -----------------------------------------------------------------
    // GET NUMBER OF NULLS PER SAMPLE
    // -----------------------------------------------------------------

    NB_MISSING_VALUES(
        ch_counts,
        ch_valid_gene_ids.collect()
    )

    ch_nb_nulls_per_sample_file = NB_MISSING_VALUES.out.nb_nulls
                                    .splitCsv( header: true )
                                    .collectFile(
                                        name: 'nb_nulls_per_sample.csv',
                                        seed: "sample,count",
                                        newLine: true,
                                        storeDir: "${outdir}/statistics/",
                                        sort: true
                                    )
                                    {
                                        item -> "${item["sample"]},${item["count"]}"
                                    }

    emit:
    nb_nulls_per_sample_file = ch_nb_nulls_per_sample_file

}
