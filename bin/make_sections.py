#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
"""
make_sections.py - Make gene sections based on a weighted average of mean gene expression across platforms
"""

import argparse
import logging
from pathlib import Path

import config
from common import write_csv_with_floats

import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
        "--nb-sections",
        type=int,
        dest="nb_sections",
        required=True,
        help="Number of sections to divide the data into",
    )
    return parser.parse_args()


def parse_stat_file(file: Path) -> pl.DataFrame:
    return pl.scan_csv(file).select(
        pl.col(config.GENE_ID_COLNAME).cast(pl.String()),
        pl.col(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME).cast(pl.Int32()),
    ).collect()


def parse_expression_level_quantiles(files: list[Path]) -> pl.DataFrame:
    dfs = [parse_stat_file(file) for file in files]
    df = dfs[0]
    if len(dfs) == 1:
        return df
    for i, other_df in enumerate(dfs[1:]):
        df = df.join(
            other_df,
            how="full",
            on=config.GENE_ID_COLNAME,
            coalesce=True,
            suffix=f"_{i+1}"
        )
    return df


def make_sections(stat_df: pl.DataFrame, nb_sections: int):
    """
    Assigns gene to sections bases on mean expression level
    Polars only ranks non-null values and preserves the null ones.
    Expression level quantiles represent the category of expression level of each gene in a certain platform.
    The higher the quantile, the higher the average expression.
    Therefore, we want to put all genes with the highest expression in section 1, and so on...
    """
    return stat_df.with_columns(
            pl.col(config.GENE_ID_COLNAME),
            pl.mean_horizontal(pl.exclude(config.GENE_ID_COLNAME)).alias("mean_expr_level_quantiles")
        ).with_columns(
            (
                pl.col("mean_expr_level_quantiles").rank(method="ordinal", descending=True)
                / pl.col("mean_expr_level_quantiles").count()
                * nb_sections + 1
            ).floor()
            .cast(pl.UInt8)
            # we want the only value at <nb_sections +1> to be at <nb_sections>
            .replace({nb_sections + 1: nb_sections})
            .alias(config.SECTION_COLNAME)
        ).select([config.GENE_ID_COLNAME, config.SECTION_COLNAME])


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing expression level quantile intervals")
    df = parse_expression_level_quantiles(args.stat_files)

    logger.info("Getting sections")
    section_df = make_sections(df, args.nb_sections)

    logger.info(f"Writing section to {OUTFILE}")
    write_csv_with_floats(section_df, OUTFILE)


if __name__ == "__main__":
    main()
