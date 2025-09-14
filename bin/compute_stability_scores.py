#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
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

    WEIGHT: ClassVar[dict] = {
        "std": 1000,
        "cv": 1000,
        "mad": 1000,
        "normfinder": 1000,
        "genorm": 1,
    }

    WEIGHT_RATIO_NB_NULLS_TO_SCORING: ClassVar[float] = 0.01

    lf: pl.LazyFrame
    scoring_base: str

    stability_base_col: str = field(init=False)
    weight_stability_base: float = field(init=False)

    def __post_init__(self):
        self.stability_base_col = config.SCORING_BASE_TO_STABILITY_SCORE_COLUMN[self.scoring_base]
        self.weight_stability_base = self.WEIGHT.get(self.scoring_base, 1)
        self.compute_stability_score()


    def compute_stability_score(self) -> pl.LazyFrame:
        logger.info("Computing stability score for candidate genes")
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # GENERAL FORMULA FOR STABILITY
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        stability_scoring_expr = (
            pl.col(self.stability_base_col) * self.weight_stability_base
            + pl.col(config.RATIO_NULLS_VALID_SAMPLES_COLNAME) * self.WEIGHT_RATIO_NB_NULLS_TO_SCORING
        )
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        expr = (
            pl.when(pl.col(config.IS_CANDIDATE_COLNAME).is_not_null())
            .then(stability_scoring_expr)
            .otherwise(None)
        )
        # add stability score column
        self.lf = self.lf.with_columns(expr.alias(config.STABILITY_SCORE_COLNAME))


    def get_statistics_with_stability_scores(self):
        return (
            self.lf
            .sort(config.STABILITY_SCORE_COLNAME, descending=False, nulls_last=True)
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
        help="Platform stat file"
    )
    parser.add_argument(
        "--stabilities",
        type=str,
        dest="stability_files",
        required=True,
        help="Output files of Normfinder / Genorm",
    )
    parser.add_argument(
        "--scoring-base",
        type=str,
        dest="scoring_base",
        required=True,
        help="Base statistical measurement (computed: Normfinder, Genorm) "
             "or descriptive (standard deviation, coefficient of variation) to use as base for stability scoring."
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
    scored_lf.collect().write_csv(STATISTICS_WITH_SCORES_OUTFILENAME)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    stability_files = [Path(file) for file in args.stability_files.split(" ")]
    stat_files = [Path(file) for file in args.platform_stat_files.split(" ")]

    # getting metadata and mappings
    stability_df = get_stabilities(stability_files)
    stat_df = get_statistics(stat_files)

    # merges base statistics with computed stability measurements
    lf = stat_df.join(stability_df, on=config.ENSEMBL_GENE_ID_COLNAME, how="left")

    # sort genes according to the metrics present in the dataframe
    stability_scorer = StabilityScorer(lf, scoring_base=args.scoring_base)
    scored_lf = stability_scorer.get_statistics_with_stability_scores()

    # exporting computed data
    export_data(scored_lf)


if __name__ == "__main__":
    main()
