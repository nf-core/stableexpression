include { GET_CANDIDATE_GENES                } from '../../../modules/local/get_candidate_genes'
include { NORMFINDER                         } from '../../../modules/local/normfinder'
include { COMPUTE_STABILITY_SCORES           } from '../../../modules/local/compute_stability_scores'

include { GENORM                             } from '../genorm'
/*
========================================================================================
    COMPUTE STABILITY SCORES
========================================================================================
*/

workflow STABILITY_SCORING {

    take:
    ch_counts
    ch_design
    ch_stats

    main:

    // -----------------------------------------------------------------
    // GETTING CANDIDATE GENES
    // -----------------------------------------------------------------

    GET_CANDIDATE_GENES(
        ch_counts,
        ch_stats,
        params.candidate_selection_descriptor,
        params.nb_top_gene_candidates
    )
    GET_CANDIDATE_GENES.out.counts.set { ch_candidate_gene_counts }

    // -----------------------------------------------------------------
    // NORMFINDER
    // -----------------------------------------------------------------

    NORMFINDER (
        ch_candidate_gene_counts,
        ch_design
    )
    NORMFINDER.out.stability_values.set { ch_stability_scores }

    // -----------------------------------------------------------------
    // GENORM
    // -----------------------------------------------------------------

    if ( !params.skip_genorm ) {
        GENORM ( ch_candidate_gene_counts )

        ch_stability_scores
            .mix ( GENORM.out.m_measures )
            .set { ch_stability_scores }
    }

    // -----------------------------------------------------------------
    // AGGREGATION AND FINAL STABILITY SCORE
    // -----------------------------------------------------------------

    COMPUTE_STABILITY_SCORES (
        ch_stats,
        ch_stability_scores.collect(),
        params.scoring_base
    )


    emit:
    summary_statistics      = COMPUTE_STABILITY_SCORES.out.stats_with_stability_scores

}
