#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import polars as pl
from common import compute_log2, export_parquet, parse_count_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


OUTFILE_SUFFIX = ".cpm.parquet"

WARNING_REASON_FILE = "warning_reason.txt"
FAILURE_REASON_FILE = "failure_reason.txt"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Normalise data to CPM")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    return parser.parse_args()


def calculate_cpm(df: pl.DataFrame) -> pl.DataFrame:
    """
    Calculate CPM (Counts Per Million) from raw count data.

    Parameters:
    -----------
    counts_df : polars.DataFrame
        DataFrame with genes as rows and samples as columns

    Returns:
    --------
    cpm_df : polars.DataFrame
        DataFrame with CPM values
    """
    # Calculate total counts per sample (column sums)
    sums = df.select(pl.exclude(config.GENE_ID_COLNAME).sum())

    # Calculate CPM: (count / total_counts) * 1,000,000
    count_columns = df.select(pl.exclude(config.GENE_ID_COLNAME)).columns
    return df.select(
        [pl.col(config.GENE_ID_COLNAME)]
        + [(pl.col(col) / sums[col][0] * 1e6).alias(col) for col in count_columns]
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing data")

    try:
        count_df = parse_count_table(args.count_file)

        logger.info(f"Normalising {args.count_file.name}")
        count_df = calculate_cpm(count_df)

        logger.info("Computing log2 values")
        count_df = compute_log2(count_df)

        export_parquet(count_df, args.count_file, OUTFILE_SUFFIX)

    except Exception as e:
        logger.error(f"Error occurred while normalising data: {e}")
        msg = "UNEXPECTED ERROR"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)


if __name__ == "__main__":
    main()
