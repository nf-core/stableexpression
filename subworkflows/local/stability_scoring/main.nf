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
    candidate_selection_descriptor
    nb_top_gene_candidates
    min_expr_threshold
    run_genorm
    stability_score_weights

    main:

    // -----------------------------------------------------------------
    // GETTING CANDIDATE GENES
    // -----------------------------------------------------------------

    GET_CANDIDATE_GENES(
        ch_counts.collect(),
        ch_stats.collect(),
        candidate_selection_descriptor,
        nb_top_gene_candidates,
        min_expr_threshold
    )
    GET_CANDIDATE_GENES.out.counts.set { ch_candidate_gene_counts }

    // -----------------------------------------------------------------
    // NORMFINDER
    // -----------------------------------------------------------------

    NORMFINDER (
        ch_candidate_gene_counts.collect(),
        ch_design.collect()
    )
    NORMFINDER.out.stability_values.set { ch_normfinder_stability }

    // -----------------------------------------------------------------
    // GENORM
    // -----------------------------------------------------------------

    if ( run_genorm ) {
        GENORM ( ch_candidate_gene_counts )
        GENORM.out.m_measures.set { ch_genorm_stability }
    } else {
        ch_genorm_stability = Channel.value([])
    }

    // -----------------------------------------------------------------
    // AGGREGATION AND FINAL STABILITY SCORE
    // -----------------------------------------------------------------

    COMPUTE_STABILITY_SCORES (
        ch_stats.collect(),
        stability_score_weights,
        ch_normfinder_stability,
        ch_genorm_stability
    )


    emit:
    summary_statistics      = COMPUTE_STABILITY_SCORES.out.stats_with_stability_scores

}
