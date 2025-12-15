//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAKE_CHUNKS                  } from '../../../modules/local/genorm/make_chunks'
include { CROSS_JOIN                   } from '../../../modules/local/genorm/cross_join'
include { EXPRESSION_RATIO             } from '../../../modules/local/genorm/expression_ratio'
include { RATIO_STANDARD_VARIATION     } from '../../../modules/local/genorm/ratio_standard_variation'
include { COMPUTE_M_MEASURE            } from '../../../modules/local/genorm/compute_m_measure'

/*
========================================================================================
    SUBWORKFLOW TO COMPUTE PAIRWISE GENE VARIATION
========================================================================================
*/

workflow GENORM {

    take:
    ch_counts


    main:

    MAKE_CHUNKS( ch_counts )

    // we need to flatten to set each chunk file as a separate item in the channel
    ch_count_chunks = MAKE_CHUNKS.out.chunks.flatten()
    getUniqueFilePairs( ch_count_chunks ) | CROSS_JOIN

    CROSS_JOIN.out.data | EXPRESSION_RATIO

    EXPRESSION_RATIO.out.data | RATIO_STANDARD_VARIATION

    COMPUTE_M_MEASURE(
        ch_counts,
        RATIO_STANDARD_VARIATION.out.data.collect( sort: true )
    )

    emit:
    m_measures = COMPUTE_M_MEASURE.out.m_measures

}



/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Generate channels consisting of unique pairs of files
// Gets
//
def getUniqueFilePairs( ch_count_chunks ) {

    def ch_count_chunks_with_indexes = ch_count_chunks
                                        .map { file -> [file.name.tokenize('.')[1], file] } // extract file index

    return ch_count_chunks_with_indexes
            .combine( ch_count_chunks_with_indexes ) // full cartesian product with itself
            .map { // steps not mandatory but helps to make the filter clearer
                index_1, file_1, index_2, file_2 ->
                    [index_1: index_1, index_2: index_2, file_1: file_1, file_2: file_2]
            }
            .filter { it -> it.index_1 <= it.index_2 } // keeps only pairs where i <= j
            .map {
                it ->
                    def meta = [index_1: it.index_1, index_2: it.index_2] // puts indexes in a meta tuple
                    [ meta, it.file_1, it.file_2 ]
            }
}
