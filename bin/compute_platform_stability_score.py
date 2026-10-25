#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

import config
import polars as pl
from common import export_parquet, zscore_to_unit

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
STATISTICS_WITH_SCORES_OUTFILENAME = "stats_with_scores.parquet"


@dataclass
class StabilityScorer:
    N_QUANTILES: ClassVar[int] = 1000

    WEIGHT_FIELDS: ClassVar[list[str]] = [
        config.NORMFINDER_STABILITY_VALUE_COLNAME,
        config.GENORM_M_MEASURE_COLNAME,
        config.COEFFICIENT_OF_VARIATION_COLNAME,
        config.ROBUST_COEFFICIENT_OF_VARIATION_MEDIAN_COLNAME,
    ]

    WEIGHT_RATIO_NB_NULLS_TO_SCORING: ClassVar[float] = 1

    df: pl.DataFrame
    stability_score_weights_str: str
    weights: dict[str, float] = field(default_factory=dict)

    def __post_init__(self):
        self.parse_stability_score_weights()
        self.compute_stability_score()

    def parse_stability_score_weights(self):
        for weight_field, weight in zip(
            self.WEIGHT_FIELDS, self.stability_score_weights_str.split(",")
        ):
            self.weights[weight_field] = float(weight)

    @staticmethod
    def get_normalised_col(col: str) -> str:
        return f"{col}_normalised"

    def compute_stability_score(self):
        logger.info("Computing stability score for candidate genes")

        # since Normfinder is always run for candidate genes
        # we can distinguish between candidate and non-candidate genes easily with this column
        self.df = self.df.with_columns(
            pl.when(pl.col(config.NORMFINDER_STABILITY_VALUE_COLNAME).is_not_null())
            .then(1)
            .otherwise(0)
            .alias(config.IS_CANDIDATE_COLNAME)
        )

        # dividing the dataframe into two parts: candidate and non-candidate genes
        candidate_df = self.df.filter(
            pl.col(config.IS_CANDIDATE_COLNAME) == 1
        )  # keep only candidate genes
        non_candidate_df = self.df.filter(pl.col(config.IS_CANDIDATE_COLNAME) == 0)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # DATA NORMALISATION (TO [0, 1])
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        normalised_data = {}
        null_data = {}
        weight_sum = 0
        # iterate over columns that can participate in stability score calculation
        cols_for_stability_scoring = []
        for col, weight in self.weights.items():

            # if a column is absent, skip it
            if col not in self.df.columns:
                logger.warning(f"Column {col} not found in dataframe")
                continue

            normalised_col = self.get_normalised_col(col)
            # creating a null column with same name
            null_data[col] = pl.Series(normalised_col, [None] * len(non_candidate_df))

            # if all values in column are nulls, skipping it
            # this may happen for instance for robust_coefficient_of_variation_median
            # in sections where almost all genes have a median at zero
            if self.df[col].is_nan().all():
                logger.warning(f"Column {col} has only null values")
                normalised_data[col] = self.df[col].alias(normalised_col)
            else:
                normalised_data[col] = zscore_to_unit(candidate_df[col]).alias(normalised_col)
                # counting the sum of weights corresponding to the columns present
                # so that we can normalise the weights afterwards
                weight_sum += weight
                # adding this column to the list of valid cols
                cols_for_stability_scoring.append(col)

        # replacing original data with normalised ones
        candidate_df = candidate_df.with_columns(
            data for data in normalised_data.values()
        )
        # adding null columns to the non-candidate df to allow concatenation
        non_candidate_df = non_candidate_df.with_columns(
            data for data in null_data.values()
        )

        # concatenating with non candidate genes to have all genes
        self.df = pl.concat([candidate_df, non_candidate_df])

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # GENERAL FORMULA FOR STABILITY
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # adding penalty for samples with null values
        stability_scoring_expr = (
            pl.col(config.RATIO_NULLS_VALID_SAMPLES_COLNAME)
            * self.WEIGHT_RATIO_NB_NULLS_TO_SCORING
        )

        for col, weight in self.weights.items():
            if col not in cols_for_stability_scoring:
                continue
            normalised_col = self.get_normalised_col(col)
            # we do not want to include null / nan values in the stability score calculation
            # because this would result in a total null / nan value for the stability score
            # since setting 0 for these samples would favour them, we add the max value of the distribution for these samples
            stability_scoring_expr += (
                pl.when(
                    pl.col(normalised_col).is_not_null()
                    & pl.col(normalised_col).is_not_nan()
                )
                .then(pl.col(normalised_col) * weight)
                .otherwise(pl.col(normalised_col).max())
            )

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        expr = (
            pl.when(pl.col(config.IS_CANDIDATE_COLNAME) == 1)
            .then(stability_scoring_expr)
            .otherwise(None)
        )
        # add stability score column
        self.df = self.df.with_columns(expr.alias(config.STABILITY_SCORE_COLNAME))

    def get_statistics_with_stability_scores(self) -> pl.DataFrame:
        return (
            self.df.sort(
                config.STABILITY_SCORE_COLNAME, descending=False, nulls_last=True
            )
            .with_row_index(name="index")
            .with_columns((pl.col("index") + 1).alias(config.RANK_COLNAME))
            .drop("index")
        )


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Computes stability score for each gene for a specific platform"
    )
    parser.add_argument(
        "--platform",
        type=str,
        required=True,
        help="Platform name",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        dest="stats_file",
        required=True,
        help="Gene Statistics file",
    )
    parser.add_argument(
        "--normfinder-stability",
        type=str,
        required=True,
        dest="normfinder_stability_file",
        help="Output files of Normfinder",
    )
    parser.add_argument(
        "--genorm-stability",
        type=str,
        dest="genorm_stability_file",
        help="Output files of Genorm",
    )
    parser.add_argument(
        "--weights",
        dest="stability_score_weights",
        type=str,
        required=True,
        help="Weights for Coefficient of Variation / Robust Coefficient of Variation on Median / Normfinder / Genorm respectively. Must be a comma-separated string. Example: 0.7,0.1,0.1,0.1",
    )
    return parser.parse_args()


def get_stabilities(stability_files: list[Path]) -> pl.DataFrame:
    """Retrieve and concatenate stability values from a list of stability files."""
    df = pl.read_csv(stability_files[0])
    if len(stability_files) > 1:
        for file in stability_files[1:]:
            new_df = pl.read_csv(file)
            df = df.join(new_df, on=config.GENE_ID_COLNAME, how="left")
    return df.with_columns(pl.lit(1).alias(config.IS_CANDIDATE_COLNAME))


def rename_columns_for_platform(scored_df: pl.DataFrame, platform: str) -> pl.DataFrame:
    return scored_df.rename({
        col: f"{col}.{platform}" for col in scored_df.columns
        if col != config.GENE_ID_COLNAME
    })


def export_data(scored_df: pl.DataFrame):
    logger.info(f"Exporting stability scores to: {STATISTICS_WITH_SCORES_OUTFILENAME}")
    export_parquet(scored_df, STATISTICS_WITH_SCORES_OUTFILENAME)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    stat_df = pl.read_parquet(args.stats_file).fill_nan(None)

    stability_files = [
        Path(file)
        for file in [args.normfinder_stability_file, args.genorm_stability_file]
        if file is not None
    ]

    # getting metadata and mappings
    stability_df = get_stabilities(stability_files)
    # merges base statistics with computed stability measurements
    df = stat_df.join(stability_df, on=config.GENE_ID_COLNAME, how="left")

    # sort genes according to the metrics present in the dataframe
    stability_scorer = StabilityScorer(df, args.stability_score_weights)
    scored_df = stability_scorer.get_statistics_with_stability_scores()

    # adding a suffix to all columns except "gene_id"
    scored_df = rename_columns_for_platform(scored_df, args.platform)

    # exporting computed data
    export_data(scored_df)


if __name__ == "__main__":
    main()
