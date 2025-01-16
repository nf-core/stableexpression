#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import polars as pl
from pathlib import Path
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"

# experimentally chosen
GENE_CHUNK_SIZE = 300

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
    parser.add_argument(
        "--task-attempts",
        dest="task_attempts",
        type=int,
        default=1,
        help="Number of task attempts",
    )
    return parser.parse_args()


def get_nb_rows(lf: pl.LazyFrame):
    return lf.select(pl.len()).collect().item()


def parse_count_dataset(file: Path, low_memory: bool) -> pl.LazyFrame:
    lf = pl.scan_parquet(file, low_memory=low_memory).fill_null(0).fill_nan(0)
    count_columns = get_count_columns(lf)
    cols = [pl.col(ENSEMBL_GENE_ID_COLNAME)] + [
        pl.col(column).replace({0: 1e-8}).cast(pl.Float32) for column in count_columns
    ]
    return lf.select(cols)


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col
        for col in lf.collect_schema().names()
        if not col.startswith(ENSEMBL_GENE_ID_COLNAME)
    ]


def split_count_summary_in_chunks(lf: pl.LazyFrame):
    lf = lf.with_row_index(name="index")
    nb_rows = get_nb_rows(lf)
    logger.info(f"Number of rows (genes) in count file: {nb_rows}")
    nb_chunks = nb_rows // GENE_CHUNK_SIZE + 1
    logger.info(f"Number of chunks: {nb_chunks}")

    for i, start in enumerate(range(0, nb_rows + 1, GENE_CHUNK_SIZE)):
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

    low_memory = True if args.task_attempts > 1 else False
    logger.info("Parsing count file")
    lf = parse_count_dataset(args.count_file, low_memory)

    logger.info("Splitting count file into chunks")
    split_count_summary_in_chunks(lf)


if __name__ == "__main__":
    main()
