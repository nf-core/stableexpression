#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path
from math import sqrt

import config
import polars as pl
from common import write_csv_with_floats

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
STATISTICS_WITH_SCORES_OUTFILENAME = "stats_with_scores.csv"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Computes global stability score for each gene"
    )
    parser.add_argument(
        "--platform-stats-scores",
        type=str,
        dest="platform_stats_scores_file",
        required=True,
        help="Platform-specific statistics / scores file",
    )
    parser.add_argument(
        "--nb-samples-per-platform",
        type=Path,
        dest="platform_size_file",
        required=True,
        help="File containing the number of samples per platform",
    )
    parser.add_argument(
        "--std-penalty-weight",
        dest="std_penalty_weight",
        type=float,
        required=True,
        help="Weight parameter for the variance penalty",
    )
    parser.add_argument(
        "--null-penalty-weight",
        dest="null_penalty_weight",
        type=float,
        required=True,
        default = 3,
        help="Weight parameter for the null values penalty",
    )
    return parser.parse_args()


def get_scores(files: list[Path]) -> pl.DataFrame:
    """Retrieve and concatenate stats and scores values from a list of files."""
    df = pl.read_csv(files[0])
    if len(files) > 1:
        for file in files[1:]:
            new_df = pl.read_csv(file)
            df = df.join(new_df, on=config.GENE_ID_COLNAME, how="full", coalesce=True)
    return df


def compute_global_score(
    df: pl.DataFrame,
    platform_sizes: dict[str, int],
    std_penalty_weight: float,
    null_penalty_weight: float
) -> pl.DataFrame:
    """
    Compute the global stability score by weighting stability scores from multiple platforms.
    """
    stability_score_columns = [col for col in df.columns if col.startswith(config.STABILITY_SCORE_COLNAME)]

    if len(stability_score_columns) == 1: # if only one platform, we just take the only score
        return df.with_columns(pl.col(stability_score_columns[0]).alias(config.GLOBAL_STABILITY_SCORE_COLNAME))

    weighted_stability_score_columns = [f"{col}.weighted" for col in stability_score_columns]
    sum_of_weights = sum([sqrt(size) for size in platform_sizes.values()])

    # 1 - multiplying each stability score column by the square root of the nb of samples in the associated platform
    # 2 - compute the average of these weighted scores
    # 3 - mitigating this sum by a penalty term, which is the std of all platform stability scores, multiplied by alpha
    return (
        df.with_columns([
            (pl.col(col) * sqrt(platform_sizes[col.split('.')[-1]])).alias(weighted_col)
            for col, weighted_col in zip(stability_score_columns, weighted_stability_score_columns)
        ]).with_columns(
            (pl.sum_horizontal(weighted_stability_score_columns) / sum_of_weights).alias('stability_score_weighted_average'),
            pl.concat_list(stability_score_columns).list.drop_nulls().list.std().fill_null(0).alias('stability_score_std'),
            pl.sum_horizontal(pl.col(stability_score_columns).is_null()).alias('nb_null_stability_scores')
        ).with_columns(
            (
                pl.col('stability_score_weighted_average')
                + std_penalty_weight * pl.col('stability_score_std')
                + null_penalty_weight * pl.col('nb_null_stability_scores')
            ).alias(config.GLOBAL_STABILITY_SCORE_COLNAME)
        )
    )


def add_global_rank(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.sort(
            config.GLOBAL_STABILITY_SCORE_COLNAME, descending=False, nulls_last=True
        )
        .with_row_index(name="index")
        .with_columns((pl.col("index") + 1).alias(config.RANK_COLNAME))
        .drop("index")
    )


def export_data(df: pl.DataFrame):
    """Export gene expression data to CSV files."""
    logger.info(f"Exporting stability scores to: {STATISTICS_WITH_SCORES_OUTFILENAME}")
    write_csv_with_floats(df, STATISTICS_WITH_SCORES_OUTFILENAME, float_precision=5)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    df = get_scores(args.platform_stats_scores_file.split(' '))

    platform_size_df = pl.read_csv(args.platform_size_file)
    platform_sizes = dict(zip(platform_size_df["platform"], platform_size_df["nb_samples"]))

    df = compute_global_score(df, platform_sizes, args.std_penalty_weight, args.null_penalty_weight)

    df = add_global_rank(df)

    # exporting computed data
    export_data(df)


if __name__ == "__main__":
    main()
