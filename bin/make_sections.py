#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
"""
make_sections.py - Make gene sections based on a weighted average of mean gene expression across platforms
"""

import argparse
import logging
from pathlib import Path
from math import sqrt

import config
from common import write_csv_with_floats

import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

RANDOM_SEED = 42

OUTFILE = "sections.csv"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get statistics from count data for each gene"
    )
    parser.add_argument(
        "--stats",
        type=Path,
        nargs="+",
        dest="stat_files",
        required=True,
        help="Files containing statistics of expression over all datasets, one for each platform",
    )
    parser.add_argument(
        "--nb-samples-per-platform",
        type=Path,
        dest="platform_size_file",
        required=True,
        help="File containing the number of samples per platform",
    )
    parser.add_argument(
        "--nb-sections",
        type=int,
        dest="nb_sections",
        required=True,
        help="Number of sections to divide the data into",
    )
    return parser.parse_args()


def parse_stat_file(file: Path, platform: str) -> pl.DataFrame:
    return pl.scan_csv(file).select(
        pl.col(config.GENE_ID_COLNAME).cast(pl.String()),
        pl.col(config.EXPRESSION_LEVEL_QUANTILE_COLNAME).cast(pl.Int32()).alias(f"{config.EXPRESSION_LEVEL_QUANTILE_COLNAME}.{platform}"),
    ).collect()


def parse_expression_level_quantiles(files: list[Path]) -> pl.DataFrame:
    """
    Parse expression level quantiles from each platform statistics file.
    for each file, the name of the platform is at the beginning of the file name.
    """
    df = None
    for file in files:
        platform = file.name.split('.')[0]
        platform_df = parse_stat_file(file, platform)
        if df is None:
            df = platform_df
        else:
            df = df.join(
                platform_df,
                how="full",
                on=config.GENE_ID_COLNAME,
                coalesce=True,
                maintain_order="left_right"
            )

    if df is None:
        raise ValueError('No stat file found!')

    return df


def make_sections(df: pl.DataFrame, platform_sizes: dict[str, int], nb_sections: int):
    """
    Assigns gene to sections bases on mean expression level
    Polars only ranks non-null values and preserves the null ones.
    Expression level quantiles represent the category of expression level of each gene in a certain platform.
    The higher the quantile, the higher the average expression.
    Therefore, we want to put all genes with the highest expression in section 1, and so on...
    """
    # storing current list of columns
    cols = [col for col in df.select(pl.exclude(config.GENE_ID_COLNAME)).columns]
    weighted_cols = [f"{col}.weighted" for col in cols]
    sum_of_weights = sum([sqrt(size) for size in platform_sizes.values()])

    ranking_kwargs = dict(
        method="random",
        seed=RANDOM_SEED,
        descending=True
    )

    return (
        df.with_columns([
            (pl.col(col) * sqrt(platform_sizes[col.split('.')[-1]])).alias(weighted_col)
            for col, weighted_col in zip(cols, weighted_cols)
        ])
        .with_columns(
            pl.col(config.GENE_ID_COLNAME),
            (pl.mean_horizontal(weighted_cols) / sum_of_weights).alias("weighted_mean_expr_level_quantiles")
        ).with_columns(
            (
                pl.col("weighted_mean_expr_level_quantiles").rank(**ranking_kwargs)
                / pl.col("weighted_mean_expr_level_quantiles").count()
                * nb_sections + 1
            ).floor()
            .cast(pl.UInt32)
            # we want the only value at <nb_sections +1> to be at <nb_sections>
            .replace({nb_sections + 1: nb_sections})
            .alias(config.SECTION_COLNAME)
        )
        .select([config.GENE_ID_COLNAME, config.SECTION_COLNAME])
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing expression level quantile intervals")
    df = parse_expression_level_quantiles(args.stat_files)

    platform_size_df = pl.read_csv(args.platform_size_file)
    platform_sizes = dict(zip(platform_size_df["platform"], platform_size_df["nb_samples"]))

    logger.info("Getting sections")
    section_df = make_sections(df, platform_sizes, args.nb_sections)

    logger.info(f"Writing section to {OUTFILE}")
    write_csv_with_floats(section_df, OUTFILE)


if __name__ == "__main__":
    main()
