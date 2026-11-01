#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import polars as pl
from common import (
    export_parquet,
    parse_count_table,
    parse_table,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


OUTFILE_SUFFIX = ".rpk.parquet"

WARNING_REASON_FILE = "warning_reason.txt"
FAILURE_REASON_FILE = "failure_reason.txt"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Normalise data to RPK")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--gene-lengths",
        type=Path,
        dest="gene_lengths_file",
        required=True,
        help="Gene lengths file (CSV format)",
    )
    return parser.parse_args()


def compute_rpk(df: pl.DataFrame, cdna_length_df: pl.DataFrame) -> pl.DataFrame:
    """
    Process raw counts to RPK.
    """
    logger.info("Computing RPK.")
    df = df.join(cdna_length_df, on=config.GENE_ID_COLNAME)
    return df.select(
        pl.col(config.GENE_ID_COLNAME),
        pl.exclude([config.GENE_ID_COLNAME, config.CDNA_LENGTH_COLNAME]).truediv(
            pl.col(config.CDNA_LENGTH_COLNAME)
        ),
    )


def parse_gene_length(file: Path) -> pl.DataFrame:
    df = parse_table(file)
    return df.with_columns(pl.col(config.CDNA_LENGTH_COLNAME).cast(pl.UInt32))


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    try:
        logger.info("Parsing data")
        count_df = parse_count_table(args.count_file)
        cdna_length_df = parse_gene_length(args.gene_lengths_file)

        # casting to Float64 to avoid inconsistency during the subsequent computations
        count_df = count_df.with_columns(
            pl.exclude(config.GENE_ID_COLNAME).cast(pl.Float64)
        )

        logger.info(f"Normalising {args.count_file.name}")
        count_df = compute_rpk(count_df, cdna_length_df)

        outfilename = args.count_file.with_suffix(OUTFILE_SUFFIX).name
        export_parquet(count_df, outfilename)

    except Exception as e:
        logger.error(f"Error occurred while normalising data: {e}")
        msg = "UNEXPECTED ERROR"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)


if __name__ == "__main__":
    main()
