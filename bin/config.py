
# general column names
ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
RANK_COLNAME = "Rank"

# base statistics
VARIATION_COEFFICIENT_COLNAME = "variation_coefficient"
STANDARD_DEVIATION_COLNAME = "standard_deviation"
STABILITY_SCORE_COLNAME = "stability_score"
MEAN_COLNAME = "mean"
MEDIAN_COLNAME = "median"
MAD_COLNAME = "median_absolute_deviation"
EXPRESSION_LEVEL_STATUS_COLNAME = "expression_level_status"
EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME = "expression_level_quantile_interval"
RATIO_NULLS_COLNAME = "ratio_nulls_in_all_samples"
RATIO_NULLS_VALID_SAMPLES_COLNAME = "ratio_nulls_in_valid_samples"
RATIO_ZEROS_COLNAME = "ratio_zeros"
IS_CANDIDATE_COLNAME = "is_candidate"

# dataset statistics
KS_TEST_COLNAME = "kolmogorov_smirnov_pvalue"

# count dataframe
GENE_COUNT_COLNAME = "count"
SAMPLE_COLNAME = "sample"

# gene metadata
ORIGINAL_GENE_ID_COLNAME = "original_gene_id"
ORIGINAL_GENE_IDS_COLNAME = "original_gene_ids"
GENE_NAME_COLNAME = "name"
GENE_DESCRIPTION_COLNAME = "description"

# computed stability values
NORMFINDER_STABILITY_VALUE_COLNAME = "normfinder_stability_value"
GENORM_M_MEASURE_COLNAME = "genorm_m_measure"
RATIOS_STD_COLNAME = "ratios_stds"

SCORING_BASE_TO_STABILITY_SCORE_COLUMN = {
    "normfinder": NORMFINDER_STABILITY_VALUE_COLNAME,
    "genorm": GENORM_M_MEASURE_COLNAME,
    "std": STANDARD_DEVIATION_COLNAME,
    "cv": VARIATION_COEFFICIENT_COLNAME,
    "mad": MAD_COLNAME
}


