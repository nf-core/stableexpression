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
    SUBWORKFLOW TO COMPUTE PAIRWISE GENE VARIATION (ADAPTED VERSION OF GENORM)
========================================================================================
*/

workflow GENORM {

    take:
    ch_counts


    main:

    // -----------------------------------------------------------------
    // MAKE CHUNKS OF GENE COUNTS
    // -----------------------------------------------------------------

    MAKE_CHUNKS( ch_counts )

    // we need to flatten to set each chunk file as a separate item in the channel
    ch_count_chunks = getUniqueFilePairs( MAKE_CHUNKS.out.chunks.transpose() )

    // -----------------------------------------------------------------
    // CROSS JOIN CHUNKS
    // -----------------------------------------------------------------

    CROSS_JOIN( ch_count_chunks )

    // -----------------------------------------------------------------
    // PAIRWISE EXPRESSION RATIOS
    // -----------------------------------------------------------------

    EXPRESSION_RATIO( CROSS_JOIN.out.data )

    // -----------------------------------------------------------------
    // STANDARD VARIATION OF EXPRESSION RATIOS
    // -----------------------------------------------------------------

    RATIO_STANDARD_VARIATION( EXPRESSION_RATIO.out.data )

    // -----------------------------------------------------------------
    // COMPUTE M-MEASURE
    // -----------------------------------------------------------------

    ch_ratio_files = RATIO_STANDARD_VARIATION.out.data
                        .map{ meta, file -> [ [ section: meta.section ], file ] }
                        .groupTuple()

    COMPUTE_M_MEASURE(
        ch_counts.join( ch_ratio_files )
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
//
def getUniqueFilePairs( ch_count_chunks ) {

    def ch_count_chunks_with_indexes = ch_count_chunks
                                        .map { meta, file -> [meta, file.name.tokenize('.')[1], file] } // extract file index

    return ch_count_chunks_with_indexes
            .combine( // full cartesian product with itself, using the meta map as key
                ch_count_chunks_with_indexes,
                by: 0
            )
            .filter {
                meta, i, file_i, j, file_j -> i <= j } // keeps only pairs where i <= j
            .map {
                meta, i, file_i, j, file_j ->
                    def new_meta = meta + [ index_1: i, index_2: j ] // puts indexes in a meta tuple
                    [ new_meta, file_i, file_j ]
            }
}
