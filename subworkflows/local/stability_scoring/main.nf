include { COMPUTE_GENE_STATISTICS            } from '../../../modules/local/compute_gene_statistics'
include { GET_CANDIDATE_GENES                } from '../../../modules/local/get_candidate_genes'
include { NORMFINDER                         } from '../../../modules/local/normfinder'
include { PLATFORM_STABILITY_SCORE           } from '../../../modules/local/stability_scores/platform'
include { GLOBAL_STABILITY_SCORE             } from '../../../modules/local/stability_scores/global'

include { GENORM                             } from '../genorm'

/*
========================================================================================
    SUBWORKFLOW TO COMPUTE STABILITY SCORES
========================================================================================
*/

workflow STABILITY_SCORING {

    take:
    ch_platform_counts_design // [ meta, count_file, design]
    ch_non_imputed_counts // [ meta, count_file]
    ch_ratio_nulls_per_sample_file
    max_null_ratio_valid_sample
    nb_candidates_per_section
    nb_sections
    skip_genorm
    stability_score_weights
    outdir

    main:

    ch_platform_counts = ch_platform_counts_design.map { meta, counts, design -> [ meta, counts ] }
    ch_platform_design = ch_platform_counts_design.map { meta, counts, design -> [ meta, design ] }

    // -----------------------------------------------------------------
    // PLATFORM-SPECIFIC STATISTICS
    // -----------------------------------------------------------------

    COMPUTE_GENE_STATISTICS(
        ch_platform_counts.join( ch_non_imputed_counts ),
        ch_ratio_nulls_per_sample_file.collect(),
        max_null_ratio_valid_sample
    )
    ch_stats = COMPUTE_GENE_STATISTICS.out.stats

    // -----------------------------------------------------------------
    // GETTING CANDIDATE GENES
    // -----------------------------------------------------------------

    GET_CANDIDATE_GENES(
        ch_platform_counts.join( ch_stats ),
        nb_candidates_per_section,
        nb_sections
    )

    ch_candidate_gene_counts = splitBySection( GET_CANDIDATE_GENES.out.section_counts )
    ch_section_stats         = splitBySection( GET_CANDIDATE_GENES.out.section_stats )

    // -----------------------------------------------------------------
    // NORMFINDER
    // -----------------------------------------------------------------

    ch_normfinder_input = ch_candidate_gene_counts.map { meta, counts -> [ meta.platform, meta, counts] }
                            .join( ch_platform_design.map { meta, design -> [ meta.platform, design ] } )
                            .map { platform, meta, counts, design -> [ meta, counts, design] }

    NORMFINDER( ch_normfinder_input )

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
    // INTERMEDIATE STABILITY SCORE FOR EACH PLATFORM
    // -----------------------------------------------------------------

    ch_platform_stability_score_input = ch_normfinder_stabilities.map { meta, file -> [ "${meta.platform}_${meta.section}", meta, file] }
                                        .join( ch_genorm_stability.map { meta, file -> [ "${meta.platform}_${meta.section}", file] } )
                                        .join( ch_section_stats.map { meta, file -> [ "${meta.platform}_${meta.section}", file] } )
                                        .map { key, meta, file1, file2, file3 -> [ meta, file1, file2, file3 ] }

    PLATFORM_STABILITY_SCORE (
        ch_platform_stability_score_input,
        stability_score_weights
    )

    ch_platform_stats_with_stability_scores = PLATFORM_STABILITY_SCORE.out.stats_with_stability_scores

    // -----------------------------------------------------------------
    // COMPUTE NUMBER OF SAMPLES IN EACH PLATFORM (FOR THE COMPUTATION OF THE WEIGHTED AVERAGE IN THE GLOBAL SCORE)
    // -----------------------------------------------------------------

    ch_nb_samples_per_platform_file = ch_platform_design.map { meta, design ->
                                                def design_content = design.splitCsv( header: true )
                                                [ platform: meta.platform, size: design_content.size() ]
                                            }.collectFile(
                                                name: 'nb_samples_per_platform.csv',
                                                seed: "platform,nb_samples",
                                                newLine: true,
                                                sort: true,
                                                storeDir: "${outdir}/design/"
                                            ) {
                                                item -> "${item.platform},${item.size}"
                                            }

    // -----------------------------------------------------------------
    // AGGREGATION AND FINAL STABILITY SCORE
    // -----------------------------------------------------------------

    ch_global_stability_score_input = ch_platform_stats_with_stability_scores
                                        .map { meta, file -> [ meta.section, meta, file] }
                                        .groupTuple()
                                        .map { section, metas, files -> [ [ section: section ], files ] }

    def lambda = 0.7
    GLOBAL_STABILITY_SCORE(
        ch_global_stability_score_input,
        ch_nb_samples_per_platform_file.first(),
        lambda
    )

    emit:
    summary_statistics = GLOBAL_STABILITY_SCORE.out.stats_with_stability_scores

}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def splitBySection( ch_files ) {
    return ch_files
            .map { meta, files ->
                // if one file, wrap it in a list
                // otherwise, the collect operator separates the file path into its components,
                def fileList = files instanceof List ? files : [files]
                fileList.collect {
                    file ->
                        [ [ platform: meta.platform, section: file.name.tokenize(".")[0] ], file ]
                    }
            }
            .flatMap{ n -> n } // turns a channel of one list of n files into a channel of n files
}
