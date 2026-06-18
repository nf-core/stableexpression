#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from math import ceil
from pathlib import Path

import config
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# experimentally chosen
GENE_CHUNK_SIZE = 100
ZERO_REPLACE_VALUE = 1e-8

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Compute M-measure for each gene")
    parser.add_argument(
        "--counts",
        type=Path,
        dest="count_file",
        required=True,
        help="File containing normalised counts for all genes and all samples",
    )
    return parser.parse_args()


def get_nb_rows(lf: pl.LazyFrame):
    return lf.select(pl.len()).collect().item()


def parse_count_dataset(file: Path) -> pl.LazyFrame:
    lf = pl.scan_parquet(file).fill_null(0).fill_nan(0)
    count_columns = get_count_columns(lf)
    cols = [pl.col(config.GENE_ID_COLNAME)] + [
        pl.col(column).replace({0: ZERO_REPLACE_VALUE}).cast(pl.Float32)
        for column in count_columns
    ]
    return lf.select(cols)


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the config.GENE_ID_COLNAME column.

    The config.GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col
        for col in lf.collect_schema().names()
        if not col.startswith(config.GENE_ID_COLNAME)
    ]


def split_count_summary_in_chunks(lf: pl.LazyFrame):
    lf = lf.with_row_index(name="index")

    nb_rows = get_nb_rows(lf)
    logger.info(f"Number of rows (genes) in count file: {nb_rows}")
    nb_chunks = ceil(nb_rows / GENE_CHUNK_SIZE)
    logger.info(f"Number of chunks: {nb_chunks}")

    for i, start in enumerate(range(0, nb_rows, GENE_CHUNK_SIZE)):
        partition = (
            lf.filter(
                (pl.col("index") >= start) & (pl.col("index") < start + GENE_CHUNK_SIZE)
            )
            .drop("index")
            .collect()
        )
        outfile = f"count_chunk.{i}.parquet"
        partition.write_parquet(outfile)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing count file")
    lf = parse_count_dataset(args.count_file)

    logger.info("Splitting count file into chunks")
    split_count_summary_in_chunks(lf)


if __name__ == "__main__":
    main()
