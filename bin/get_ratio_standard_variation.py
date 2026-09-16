#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
from common import get_nb_rows

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


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the config.GENE_ID_COLNAME column.

    The config.GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col
        for col in lf.collect_schema().names()
        if not col.startswith(config.GENE_ID_COLNAME)
    ]


def compute_standard_deviations(file: Path) -> pl.LazyFrame:
    ratios_lf = pl.scan_parquet(file)
    ratio_columns = [
        col for col in ratios_lf.collect_schema().names() if col.endswith("_log_ratio")
    ]
    concat_ratios_lf = ratios_lf.select(
        [
            pl.concat_list(
                [pl.col(col) for col in ratio_columns[i : i + RATIO_CHUNK_SIZE]]
            ).alias(f"concat_list_chunk_{i // RATIO_CHUNK_SIZE}")
            for i in range(0, len(ratio_columns), RATIO_CHUNK_SIZE)
        ]
    ).select(pl.concat_list(pl.all()).alias("ratios"))
    return pl.concat(
        [
            concat_ratios_lf.select("ratios"),
            ratios_lf.select(pl.exclude("^.*_log_ratio$")),  # gene_id & gene_id_other
        ],
        how="horizontal",
    ).select(
        pl.col("ratios").list.std(ddof=0).alias(config.RATIOS_STD_COLNAME),
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

    logger.info(f"Computing standard deviations for {str(args.ratio_file)}")
    std_lf = compute_standard_deviations(args.ratio_file)
    std_lf = group_standard_deviations(std_lf)

    # when the ratio file corresponds to the same gene ids cross joined with themselves (i == i)
    # then we want only only one row per gene id
    if get_nb_rows(std_lf) == 0:
        raise ValueError(
            f"No output following treatment of file {str(args.ratio_file)}"
        )

    # sort items in each list
    std_lf = std_lf.with_columns(pl.col(config.RATIOS_STD_COLNAME).list.sort())

    outfile = args.ratio_file.name.replace("ratios", "std")
    std_lf.sink_parquet(outfile)

    logger.info(f"Wrote standard deviations to {outfile}")


if __name__ == "__main__":
    main()
