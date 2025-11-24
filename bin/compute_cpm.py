#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


CPM_NORM_SUFFIX = ".cpm.csv"

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


def parse_counts(file: Path):
    if file.suffix == ".csv":
        return pd.read_csv(file, header=0, index_col=0)
    else:  # .tsv
        return pd.read_csv(file, header=0, sep="\t", index_col=0)


def calculate_cpm(counts_df: pd.DataFrame):
    """
    Calculate CPM (Counts Per Million) from raw count data.

    Parameters:
    -----------
    counts_df : pandas.DataFrame
        DataFrame with genes as rows and samples as columns

    Returns:
    --------
    cpm_df : pandas.DataFrame
        DataFrame with CPM values
    """
    # Calculate total counts per sample (column sums)
    total_counts = counts_df.sum(axis=0)

    # Calculate CPM: (count / total_counts) * 1,000,000
    cpm_df = (counts_df / total_counts) * 1e6

    return cpm_df


def export_normalised_data(count_df: pd.DataFrame, count_file: Path):
    """Export gene expression data to CSV."""
    # replace .csv / .tsv by .tpm.csv
    outfilename = ".".join(count_file.name.split(".")[:-1]) + CPM_NORM_SUFFIX
    logger.info(f"Exporting CPM normalised counts to: {outfilename}")
    count_df.to_csv(outfilename, index=True, header=True)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing data")

    try:
        count_df = parse_counts(args.count_file)
        count_df.index.name = config.GENE_ID_COLNAME

        logger.info(f"Normalising {args.count_file.name}")

        count_df = calculate_cpm(count_df)

        export_normalised_data(count_df, args.count_file)

    except Exception as e:
        logger.error(f"Error occurred while normalising data: {e}")
        msg = "UNEXPECTED ERROR"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)


if __name__ == "__main__":
    main()
