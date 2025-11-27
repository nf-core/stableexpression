#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

COL_TO_OUTFILE = {"skewness": "skewness.txt", "ratio_zeros": "ratio_zeros.txt"}


# ALLOWED_TARGET_DISTRIBUTIONS = ["normal", "uniform"]


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute general statistics from count data for each sample"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    return parser.parse_args()


def compute_dataset_statistics(count_df: pd.DataFrame) -> pd.DataFrame:
    skewness = count_df.skew()
    ratio_zeros = (count_df == 0).sum() / len(count_df)
    return pd.DataFrame({"skewness": skewness, "ratio_zeros": ratio_zeros}).T


def export_count_data(dataset_stats_df: pd.DataFrame):
    """
    Export dataset statistics to CSV files.
    Write each statistic to a separate file, on a single row
    """
    for col, outfile_name in COL_TO_OUTFILE.items():
        logger.info(f"Exporting dataset statistics {col} to: {outfile_name}")
        pd.DataFrame(dataset_stats_df.loc[col]).T.to_csv(
            outfile_name, index=False, header=False, float_format="%.4f"
        )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_file = args.count_file

    logger.info(f"Computing dataset statistics for {count_file.name}")
    count_df = pd.read_csv(count_file, index_col=0, header=0)

    dataset_stats_df = compute_dataset_statistics(count_df)

    export_count_data(dataset_stats_df)


if __name__ == "__main__":
    main()
