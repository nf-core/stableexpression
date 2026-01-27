include { GET_NB_NULLS_PER_SAMPLE as GET_NB_NULLS              } from '../../../modules/local/get_nb_nulls_per_sample'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow TAG_SAMPLES {

    take:
    ch_all_counts      // [ [ platform: platform, dataset_size: size], file ]

    main:

    // -----------------------------------------------------------------
    // PLATFORM-SPECIFIC STATISTICS
    // -----------------------------------------------------------------

    GET_NB_NULLS( ch_all_counts.collect() )

    emit:
    nb_nulls             = GET_NB_NULLS.out.nb_nulls

}
