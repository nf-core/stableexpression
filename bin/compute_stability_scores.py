#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
from sklearn.preprocessing import QuantileTransformer
from dataclasses import dataclass, field
from typing import ClassVar
import logging

import config

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
STATISTICS_WITH_SCORES_OUTFILENAME = "stats_with_scores.csv"


@dataclass
class StabilityScorer:
    N_QUANTILES: ClassVar[int] = 1000

    WEIGHT_FIELDS: ClassVar[list] = [
        config.VARIATION_COEFFICIENT_COLNAME,
        config.MAD_COLNAME,
        config.NORMFINDER_STABILITY_VALUE_COLNAME,
        config.GENORM_M_MEASURE_COLNAME,
    ]

    WEIGHT_RATIO_NB_NULLS_TO_SCORING: ClassVar[float] = 1

    df: pl.DataFrame
    stability_score_weights_str: str
    weights: dict[str, float] = field(default_factory=dict)

    def __post_init__(self):
        self.compute_stability_score()
        self.parse_stability_score_weights()

    def parse_stability_score_weights(self):
        for field, weight in zip(
            self.WEIGHT_FIELDS, self.stability_score_weights_str.split(",")
        ):
            self.weights[field] = float(weight)

    @staticmethod
    def quantile_normalise(data: pl.Series, new_name: str) -> pl.Series:
        """
        Quantile normalize a series
        """
        array = data.to_numpy().reshape(-1, 1)
        transformer = QuantileTransformer(output_distribution="uniform")
        normalised_array = transformer.fit_transform(array)
        return pl.Series(new_name, normalised_array.ravel())

    def compute_stability_score(self) -> pl.LazyFrame:
        logger.info("Computing stability score for candidate genes")

        candidate_df = self.df.filter(
            pl.col(config.IS_CANDIDATE_COLNAME) == 1
        )  # keep only candidate genes
        non_candidate_df = self.df.filter(pl.col(config.IS_CANDIDATE_COLNAME).is_null())

        normalised_data = {}
        null_data = {}
        weight_sum = 0
        for col, weight in self.weights.items():
            if col not in self.df.columns:
                continue
            data = candidate_df.select(col).to_series()
            normalised_col = f"{col}_normalised"
            normalised_data[col] = self.quantile_normalise(
                data, new_name=normalised_col
            )
            # creating a null column with same name
            null_data[col] = pl.Series(normalised_col, [None] * len(non_candidate_df))
            # if this column is present, add its weight to the sum
            weight_sum += weight

        # replacing original data with quantile normalised ones
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
        # genes with at least one zero value are already excluded at that stage
        stability_scoring_expr = (
            pl.col(config.RATIO_NULLS_VALID_SAMPLES_COLNAME)
            * self.WEIGHT_RATIO_NB_NULLS_TO_SCORING
        )
        for col, weight in self.weights.items():
            if col not in self.df.columns:
                continue
            normalised_col = f"{col}_normalised"
            stability_scoring_expr += pl.col(normalised_col) * weight / weight_sum

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        expr = (
            pl.when(pl.col(config.IS_CANDIDATE_COLNAME).is_not_null())
            .then(stability_scoring_expr)
            .otherwise(None)
        )
        # add stability score column
        self.df = self.df.with_columns(expr.alias(config.STABILITY_SCORE_COLNAME))
        print(self.df)

    def get_statistics_with_stability_scores(self):
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
        description="Computes stability score for each gene"
    )
    parser.add_argument(
        "--stats",
        type=str,
        dest="platform_stat_files",
        required=True,
        help="Platform stat file",
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
        help="Weights for Standard deviation / Median absolute deviation / Normfinder / Genorm respectively. Must be a comma-separated string. Example: 0.7,0.1,0.1,0.1",
    )
    return parser.parse_args()


def get_stabilities(stability_files: list[Path]) -> pl.LazyFrame:
    """Retrieve and concatenate stability values from a list of stability files."""
    lf = pl.scan_csv(stability_files[0])
    if len(stability_files) > 1:
        for file in stability_files[1:]:
            new_df = pl.scan_csv(file)
            lf = lf.join(new_df, on=config.ENSEMBL_GENE_ID_COLNAME, how="left")
    return lf.with_columns(pl.lit(1).alias(config.IS_CANDIDATE_COLNAME))


def get_statistics(stat_files: list[Path]) -> pl.LazyFrame:
    """Retrieve and concatenate data from a list of statistics files."""
    lf = pl.scan_csv(stat_files[0])
    if len(stat_files) > 1:
        for file in stat_files[1:]:
            new_df = pl.scan_csv(file)
            lf = lf.join(new_df, on=config.ENSEMBL_GENE_ID_COLNAME, how="left")
    return lf


def export_data(scored_lf: pl.LazyFrame):
    """Export gene expression data to CSV files."""
    logger.info(f"Exporting stability scores to: {STATISTICS_WITH_SCORES_OUTFILENAME}")
    scored_lf.write_csv(STATISTICS_WITH_SCORES_OUTFILENAME)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    stat_files = [Path(file) for file in args.platform_stat_files.split(" ")]
    stat_lf = get_statistics(stat_files)

    stability_files = [
        Path(file)
        for file in [args.normfinder_stability_file, args.genorm_stability_file]
        if file is not None
    ]

    # getting metadata and mappings
    stability_lf = get_stabilities(stability_files)
    # merges base statistics with computed stability measurements
    lf = stat_lf.join(stability_lf, on=config.ENSEMBL_GENE_ID_COLNAME, how="left")

    # sort genes according to the metrics present in the dataframe
    stability_scorer = StabilityScorer(lf.collect(), args.stability_score_weights)
    scored_lf = stability_scorer.get_statistics_with_stability_scores()

    # exporting computed data
    export_data(scored_lf)


if __name__ == "__main__":
    main()
