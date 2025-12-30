#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".filtered.csv"

MAX_RATIO_ZEROS = 0.9


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
    return parser.parse_args()


def parse_counts(file: Path):
    if file.suffix == ".csv":
        return pd.read_csv(file, header=0, index_col=0)
    else:  # .tsv
        return pd.read_csv(file, header=0, sep="\t", index_col=0)


def filter_out_columns_with_high_zero_ratio(df: pd.DataFrame, max_ratio_zeros: float):
    zero_ratio = df.eq(0).mean(axis=0)
    return df.loc[:, zero_ratio <= max_ratio_zeros]


def export_data(df: pd.DataFrame, outfile: Path):
    logger.info(f"Exporting filtered counts to: {outfile}")
    df.to_csv(outfile, index=True, header=True)
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
    count_df = parse_counts(args.count_file)
    logger.info(
        f"Loaded count data with {len(count_df)} rows and {count_df.shape[1]} columns"
    )

    valid_count_df = filter_out_columns_with_high_zero_ratio(count_df, MAX_RATIO_ZEROS)
    if valid_count_df.shape[1] == 0:
        logger.error("No valid columns remaining")
        sys.exit(0)
    else:
        logger.info(
            f"Filtered out {count_df.shape[1] - valid_count_df.shape[1]} columns"
        )
        outfile = args.count_file.with_suffix(OUTFILE_SUFFIX)
        export_data(count_df, outfile)


if __name__ == "__main__":
    main()
