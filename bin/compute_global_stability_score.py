#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path
from math import sqrt

import config
import polars as pl
from common import write_csv_with_floats, normalise_with_z_score

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
    null_penalty_weight: float,
) -> pl.DataFrame:
    """
    Compute the global stability score by weighting stability scores from multiple platforms.

    Steps:
      1. Z-score-normalise each platform's stability score (shifted so the min is 0).
      2. Weight each score by sqrt(platform sample size) and average them.
      3. Subtract penalties: cross-platform score dispersion and missing platforms.
    """
    original_columns = df.columns
    score_cols = [
        col for col in df.columns
        if col.startswith(config.STABILITY_SCORE_COLNAME)
    ]

    # With a single platform, its score is the global score as-is.
    if len(score_cols) == 1:
        return df.with_columns(
            pl.col(score_cols[0]).alias(config.GLOBAL_STABILITY_SCORE_COLNAME)
        )

    # --- Step 1: normalise scores so platforms are comparable ---
    # Z-score per platform, then subtract the min so every platform starts at 0.
    normalised_cols = []
    for col in score_cols:
        normalised_col = f"{col}.normalised"
        df = df.hstack(
            [normalise_with_z_score(df[col], new_name=normalised_col)]
        )
        normalised_cols.append(normalised_col)

    # --- Step 2: weight by platform size and average ---
    # after normalising, col names are like 'a_string.<platform>.normalised
    # so taking col.split(".")[-2] gives the associated platform
    weights = {
        col: sqrt(platform_sizes[col.split(".")[-2]]) for col in normalised_cols
    }
    sum_of_weights = sum(weights.values())

    weighted_cols = [
        (pl.col(normalised_col) * weight).alias(f"{normalised_col}.weighted")
        for normalised_col, weight in zip(normalised_cols, weights.values())
    ]

    # --- Step 3: combine and apply penalties ---
    weighted_average = (
        pl.sum_horizontal(f"{col}.weighted" for col in normalised_cols)
        / sum_of_weights
    ).alias("stability_score_weighted_average")

    score_dispersion = (
        pl.concat_list(score_cols)
        .list.drop_nulls()
        .list.std()
        .fill_null(0)
        .alias("stability_score_std")
    )

    missing_platforms = (
        pl.sum_horizontal(pl.col(score_cols).is_null())
        .alias("nb_null_stability_scores")
    )

    global_score = (
        pl.col("stability_score_weighted_average")
        + std_penalty_weight * pl.col("stability_score_std")
        + null_penalty_weight * pl.col("nb_null_stability_scores")
    )

    return (
        df.with_columns(weighted_cols)
        .with_columns(weighted_average, score_dispersion, missing_platforms)
        .with_columns(global_score.alias(config.GLOBAL_STABILITY_SCORE_COLNAME))
        .select(original_columns + [config.GLOBAL_STABILITY_SCORE_COLNAME])
    )


def add_global_rank(df: pl.DataFrame) -> pl.DataFrame:
    """

    """
    return (
        df.sort(
            [config.GLOBAL_STABILITY_SCORE_COLNAME, config.GENE_ID_COLNAME],
            descending=[False, False],
            nulls_last=True
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
