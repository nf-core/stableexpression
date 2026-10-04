/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GET_PUBLIC_ACCESSIONS                  } from '../subworkflows/local/get_public_accessions'
include { DOWNLOAD_PUBLIC_DATASETS               } from '../subworkflows/local/download_public_datasets'
include { ID_MAPPING                             } from '../subworkflows/local/idmapping'
include { SAMPLE_FILTERING                       } from '../subworkflows/local/sample_filtering'
include { NORMALISATION                          } from '../subworkflows/local/normalisation'
include { DATASET_ANALYSIS                       } from '../subworkflows/local/dataset_analysis'
include { GENE_STATISTICS                        } from '../subworkflows/local/gene_statistics'
include { STABILITY_SCORING                      } from '../subworkflows/local/stability_scoring'
//include { REPORTING                              } from '../subworkflows/local/reporting'

include { checkCounts                            } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STABLEEXPRESSION {

    take:
    ch_input_datasets


    main:

    ch_accessions                          = channel.empty()
    ch_downloaded_datasets                 = channel.empty()
    ch_counts_ids_filtered_renamed         = channel.empty()
    ch_counts_samples_filtered             = channel.empty()
    ch_counts_first_normalissation         = channel.empty()
    ch_normalised_counts                   = channel.empty()
    ch_gene_length_file                    = channel.empty()
    ch_all_counts                          = channel.empty()
    ch_imputed_counts                      = channel.empty()
    ch_stats_all_genes_with_scores         = channel.empty()
    ch_platform_statistics                 = channel.empty()
    ch_whole_gene_metadata                 = channel.empty()
    ch_whole_gene_id_mapping               = channel.empty()

    def species = params.species.split(' ').join('_').toLowerCase()

    // -----------------------------------------------------------------
    // FETCH PUBLIC ACCESSIONS
    // -----------------------------------------------------------------

    GET_PUBLIC_ACCESSIONS(
        species,
        params.skip_fetch_eatlas_accessions,
        params.fetch_geo_accessions,
        params.platform,
        params.keywords,
        params.accessions ? channel.fromList( params.accessions.tokenize(',') ) : channel.empty(),
        params.accessions_file ? channel.fromPath(params.accessions_file, checkIfExists: true) : channel.empty(),
        params.excluded_accessions ? channel.fromList( params.excluded_accessions.tokenize(',') ) : channel.empty(),
        params.excluded_accessions_file ? channel.fromPath(params.excluded_accessions_file, checkIfExists: true) : channel.empty(),
        params.random_sampling_size,
        params.random_sampling_seed,
        params.outdir
        )

    ch_accessions = GET_PUBLIC_ACCESSIONS.out.accessions

    // -----------------------------------------------------------------
    // DOWNLOAD GEO DATASETS IF NEEDED
    // -----------------------------------------------------------------

    if ( !params.accessions_only) {

        DOWNLOAD_PUBLIC_DATASETS (
            species,
            ch_accessions
        )

        ch_downloaded_datasets = DOWNLOAD_PUBLIC_DATASETS.out.datasets

    }

    if ( !params.accessions_only && !params.download_only ) {

        ch_counts = ch_input_datasets.mix( ch_downloaded_datasets )
        // returns an error with a message if no dataset was found
        checkCounts( ch_counts, params.fetch_geo_accessions )

        // -----------------------------------------------------------------
        // IDMAPPING
        // -----------------------------------------------------------------

        // tries to map gene IDs to Ensembl IDs whenever possible
        ID_MAPPING(
            ch_counts,
            species,
            params.skip_id_mapping,
            params.skip_cleaning_gene_ids,
            params.gprofiler_target_db,
            params.gene_id_mapping,
            params.gene_metadata,
            params.min_occurrence_freq,
            params.min_occurrence_quantile,
            params.outdir
        )

        ch_counts_ids_filtered_renamed    = ID_MAPPING.out.counts
        ch_whole_gene_id_mapping          = ID_MAPPING.out.mapping
        ch_whole_gene_metadata            = ID_MAPPING.out.metadata
        ch_valid_gene_ids                 = ID_MAPPING.out.valid_gene_ids

        ch_counts = ch_counts_ids_filtered_renamed

        // -----------------------------------------------------------------
        // FILTER OUT SAMPLES NOT VALID
        // -----------------------------------------------------------------

        SAMPLE_FILTERING (
            ch_counts,
            ch_valid_gene_ids,
            params.max_zero_ratio,
            params.max_null_ratio,
            params.outdir
        )

        ch_counts_samples_filtered     = SAMPLE_FILTERING.out.counts
        ch_ratio_nulls_per_sample_file = SAMPLE_FILTERING.out.ratio_nulls_per_sample_file

        // -----------------------------------------------------------------
        // ANALYSIS OF NORMALISED DATASETS
        // -----------------------------------------------------------------

        DATASET_ANALYSIS(
            ch_counts_samples_filtered
        )

        // -----------------------------------------------------------------
        // NORMALISATION OF COUNTS
        // -----------------------------------------------------------------

        NORMALISATION(
            ch_counts_samples_filtered,
            species,
            ch_valid_gene_ids,
            params.skip_gene_length_normalisation,
            params.missing_value_imputer,
            params.quantile_normalisation,
            params.quantile_norm_target_distrib,
            params.gff,
            params.gff_url,
            params.gene_length,
            params.outdir
        )

        ch_normalised_counts          = NORMALISATION.out.normalised
        ch_imputed_counts             = NORMALISATION.out.imputed
        ch_non_imputed_counts         = NORMALISATION.out.non_imputed
        ch_gene_length_file           = NORMALISATION.out.gene_length_file

        // -----------------------------------------------------------------
        // COMPUTE BASE STATISTICS FOR ALL GENES,
        // GET CANDIDATES AS REFERENCE GENE AND COMPUTES VARIOUS STABILITY VALUES
        // -----------------------------------------------------------------

        STABILITY_SCORING (
            ch_normalised_counts,
            ch_non_imputed_counts,
            ch_ratio_nulls_per_sample_file,
            params.max_null_ratio_valid_sample,
            params.nb_candidates_per_section,
            params.nb_sections,
            params.skip_genorm,
            params.stability_score_weights
        )

        ch_stats_all_genes_with_scores = STABILITY_SCORING.out.summary_statistics

    }

    // -----------------------------------------------------------------
    // REPORTING
    // -----------------------------------------------------------------
/*
    REPORTING(
        ch_normalised_counts,
        ch_stats_all_genes_with_scores,
        ch_whole_gene_metadata,
        ch_whole_gene_id_mapping,
        params.target_genes,
        params.target_gene_file,
        params.skip_dash_app,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
    )
*/
    emit:
    accessions                             = GET_PUBLIC_ACCESSIONS.out.raw_accessions
    downloaded                             = ch_downloaded_datasets
    id_filtered_renamed                    = ch_counts_ids_filtered_renamed
    samples_filtered                       = ch_counts_samples_filtered
    normalised                             = ch_normalised_counts
    gene_length_file                       = ch_gene_length_file
    imputed                                = ch_imputed_counts
    //all_genes_summary                      = REPORTING.out.all_genes_summary
    //multiqc_report                         = REPORTING.out.multiqc_report.toList()
    //dash_app                               = REPORTING.out.dash_app
all_genes_summary = channel.empty()
multiqc_report = channel.empty().toList()
dash_app = channel.empty()
}
