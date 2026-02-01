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
    nb_candidates_per_section
    nb_sections
    skip_genorm
    stability_score_weights

    main:

    // -----------------------------------------------------------------
    // GETTING CANDIDATE GENES
    // -----------------------------------------------------------------

    GET_CANDIDATE_GENES(
        ch_counts.collect(), // single item
        ch_stats.collect(), // single item
        nb_candidates_per_section,
        nb_sections
    )
    ch_candidate_gene_counts = GET_CANDIDATE_GENES.out.counts
                                .map{ files -> ["key", files] }
                                .transpose()
                                .map {
                                    key, file ->
                                            section = file.name.tokenize(".")[0]
                                            [ [ section: section ], file]
                                }

    // -----------------------------------------------------------------
    // NORMFINDER
    // -----------------------------------------------------------------

    NORMFINDER (
        ch_candidate_gene_counts,
        ch_design.collect() // single item
    )
    ch_normfinder_stabilities = NORMFINDER.out.stability_values

    // -----------------------------------------------------------------
    // GENORM
    // -----------------------------------------------------------------

    if ( !skip_genorm ) {
        GENORM ( ch_candidate_gene_counts )
        ch_genorm_stability = GENORM.out.m_measures
    } else {
        ch_genorm_stability = channel.value([:])
    }

    // -----------------------------------------------------------------
    // AGGREGATION AND FINAL STABILITY SCORE
    // -----------------------------------------------------------------

    COMPUTE_STABILITY_SCORES (
        ch_normfinder_stabilities.join(ch_genorm_stability),
        ch_stats.collect(), // single item
        stability_score_weights
    )

    emit:
    summary_statistics      = COMPUTE_STABILITY_SCORES.out.stats_with_stability_scores

}
