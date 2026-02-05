#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import parse_count_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

KEY_TO_OUTFILE = {"skewness": "skewness.txt"}


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


def compute_dataset_statistics(df: pl.DataFrame) -> dict:
    # sample count skewness
    skewness = df.select(pl.exclude(config.GENE_ID_COLNAME).skew()).row(0)
    return dict(skewness=list(skewness))


def export_count_data(stats: dict):
    """
    Export dataset statistics to CSV files.
    Write each statistic to a separate file, on a single row
    """
    for key, outfile_name in KEY_TO_OUTFILE.items():
        logger.info(f"Exporting dataset statistics {key} to: {outfile_name}")
        with open(outfile_name, "w") as outfile:
            outfile.write(",".join([str(val) for val in stats[key]]))


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    count_file = args.count_file

    logger.info(f"Computing dataset statistics for {count_file.name}")
    count_df = parse_count_table(count_file)

    stat_dict = compute_dataset_statistics(count_df)

    export_count_data(stat_dict)


if __name__ == "__main__":
    main()
