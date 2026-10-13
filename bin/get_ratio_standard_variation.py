#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
"""
get_ratio_standard_variation.py

The input is like:
    gene_id	gene_id_other	sampleA	sampleB
    AT1G01060	AT4G34290	0.48072526	-0.049063634
    AT1G01060	AT4G34460	0.21205564	0.15575644

The output is like:
    gene_id ratio_stds
    AT1G01060   [0.04122424, 0.049532797, 0.050319236, ...]
"""

import argparse
import logging
from pathlib import Path

import config
from common import get_nb_rows, get_count_columns

import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# experimentally chosen
RATIO_CHUNK_SIZE = 100


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Compute M-measure for each gene")
    parser.add_argument(
        "--file",
        type=Path,
        dest="ratio_file",
        required=True,
        help="File log of pairwise expression ratios",
    )
    return parser.parse_args()


def compute_standard_deviations(ratios_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    For each gene pair, combine all the *_log_ratio columns into one list and computes the standard deviation across samples.
    The result is one row per pair: (ratios_std, gene_id, gene_id_other)
    """
    ratio_columns = [
        c for c in ratios_lf.collect_schema().names() if c.endswith("_log_ratio")
    ]
    if not ratio_columns:
        raise ValueError("No '*_log_ratio' columns found in input")

    # Concatenate in chunks to keep the expression tree small
    # when there are many samples
    chunks = [
        pl.concat_list(ratio_columns[i : i + RATIO_CHUNK_SIZE])
        for i in range(0, len(ratio_columns), RATIO_CHUNK_SIZE)
    ]

    return ratios_lf.select(
        pl.concat_list(chunks).list.std(ddof=0).alias(config.RATIOS_STD_COLNAME),
        pl.col(config.GENE_ID_COLNAME),
        pl.col(f"{config.GENE_ID_COLNAME}_other"),
    )


def get_column_standard_deviations(std_lf: pl.LazyFrame, column: str) -> pl.LazyFrame:
    # column is either config.GENE_ID_COLNAME or f"{config.GENE_ID_COLNAME}_other"
    return (
        std_lf.group_by(column)
        .agg(config.RATIOS_STD_COLNAME)  # getting list of ratio std for this gene
        .select(
            pl.col(column).alias(config.GENE_ID_COLNAME),
            pl.col(config.RATIOS_STD_COLNAME),
        )
    )


def group_standard_deviations(std_lf: pl.LazyFrame) -> pl.LazyFrame:
    # getting the standard devs for genes in the gene_id column
    std_a = get_column_standard_deviations(std_lf, column=config.GENE_ID_COLNAME)
    # getting the standard devs for genes in the gene_id_other column
    std_b = get_column_standard_deviations(
        std_lf, column=f"{config.GENE_ID_COLNAME}_other"
    )
    # concatenating both dataframes vertically
    # if both lists of gene ids are the identical,
    # we need to collect values only for one column to avoid duplicates
    return (
        pl.concat([std_a, std_b], how="vertical")
        .unique(subset=config.GENE_ID_COLNAME)
        .sort(
            config.GENE_ID_COLNAME
        )  # only needed to have consistent output (for snapshots)
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    ratios_lf = pl.scan_parquet(args.ratio_file)

    logger.info(f"Computing standard deviations for {args.ratio_file}")
    std_lf = compute_standard_deviations(ratios_lf)
    std_lf = group_standard_deviations(std_lf)

    # when the ratio file corresponds to the same gene ids cross joined with themselves (i == i)
    # then we want only only one row per gene id
    if get_nb_rows(std_lf) == 0:
        raise ValueError(
            f"No output following treatment of file {args.ratio_file}"
        )

    # sort items in each list
    std_lf = std_lf.with_columns(pl.col(config.RATIOS_STD_COLNAME).list.sort())

    outfile = args.ratio_file.name.replace("ratios", "std")
    std_lf.sink_parquet(outfile)

    logger.info(f"Wrote standard deviations to {outfile}")


if __name__ == "__main__":
    main()
