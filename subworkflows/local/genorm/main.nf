//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    //  CROSS JOIN SECTION DATA WITH THEMSELVES AND PAIRWISE EXPRESSION RATIOS
    // -----------------------------------------------------------------

    EXPRESSION_RATIO( ch_counts )

    // -----------------------------------------------------------------
    // STANDARD VARIATION OF EXPRESSION RATIOS
    // -----------------------------------------------------------------

    RATIO_STANDARD_VARIATION( EXPRESSION_RATIO.out.data )

    // -----------------------------------------------------------------
    // COMPUTE M-MEASURE
    // -----------------------------------------------------------------

    COMPUTE_M_MEASURE( RATIO_STANDARD_VARIATION.out.data )

    emit:
    m_measures = COMPUTE_M_MEASURE.out.m_measures

}
