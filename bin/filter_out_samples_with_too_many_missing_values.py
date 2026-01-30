#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import polars as pl
from common import export_parquet, parse_count_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".nulls_filtered.parquet"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Filter out samples not valid")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--max-null-ratio",
        type=float,
        dest="max_null_ratio",
        required=True,
        help="Maximum ratio of null values",
    )
    return parser.parse_args()


def filter_out_columns_with_high_missing_values_ratio(df: pl.DataFrame, max_null_ratio: float):
    null_ratio_df = df.select(pl.exclude(config.GENE_ID_COLNAME).is_null()).mean()
    valid_null_ratio_samples = [
        col for col in null_ratio_df.columns if null_ratio_df[col][0] <= max_null_ratio
    ]
    return df.select(pl.col(config.GENE_ID_COLNAME), pl.col(valid_null_ratio_samples))


def export_data(df: pl.DataFrame, outfile: Path):
    logger.info(f"Exporting filtered counts to: {outfile}")
    df.write_parquet(outfile)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    # putting all counts into a single dataframe
    logger.info("Loading count data...")
    count_df = parse_count_table(args.count_file)
    logger.info(
        f"Loaded count data with {len(count_df)} rows and {count_df.shape[1]} columns"
    )

    valid_count_df = filter_out_columns_with_high_missing_values_ratio(count_df, args.max_null_ratio)

    if valid_count_df.shape[1] == 0:
        logger.error("No valid columns remaining")
        sys.exit(0)
    else:
        logger.info(
            f"Filtered out {count_df.shape[1] - valid_count_df.shape[1]} columns"
        )
        export_parquet(valid_count_df, args.count_file, OUTFILE_SUFFIX)


if __name__ == "__main__":
    main()
