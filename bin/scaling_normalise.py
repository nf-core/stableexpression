#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import export_parquet, parse_count_table
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".scaled.parquet"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Z-score normalise count data for each sample in the dataset"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    return parser.parse_args()


def zscore_normalise(df: pl.DataFrame):
    """
    Z-score normalize a dataframe; column by column (each sample independently).
    """
    scaler = StandardScaler()
    return df.with_columns(
        pl.exclude(config.GENE_ID_COLNAME).map_batches(
            lambda x: scaler.fit_transform(x.to_frame()).flatten(),
            return_dtype=pl.Float32,
        )
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info(f"Parsing {args.count_file.name}")
    count_df = parse_count_table(args.count_file)

    logger.info(f"Z-score normalising {args.count_file.name}")
    zscore_normalized_counts = zscore_normalise(count_df)

    outfilename = args.count_file.with_suffix(OUTFILE_SUFFIX).name
    export_parquet(zscore_normalized_counts, outfilename)


if __name__ == "__main__":
    main()
