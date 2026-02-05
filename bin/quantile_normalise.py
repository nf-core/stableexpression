#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import export_parquet, parse_count_table
from resource_management import set_max_resources
from sklearn.preprocessing import quantile_transform

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".quant_norm.parquet"

N_QUANTILES = 1000

ALLOWED_TARGET_DISTRIBUTIONS = ["normal", "uniform"]


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Quantile normalise count data for each sample in the dataset"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--target-distrib",
        type=str,
        dest="target_distribution",
        required=True,
        choices=ALLOWED_TARGET_DISTRIBUTIONS,
        help="Target distribution to map counts to",
    )
    parser.add_argument(
        "--cpus", type=int, dest="nb_cpus", required=True, help="Number of CPUs"
    )
    parser.add_argument(
        "--memory", type=str, dest="memory", required=True, help="Memory in GB"
    )
    return parser.parse_args()


def quantile_normalise(df: pl.DataFrame, target_distribution: str):
    """
    Quantile normalize a dataframe; column by column, based on a target distribution.
    """
    kwargs = dict(
        n_quantiles=N_QUANTILES, output_distribution=target_distribution, subsample=None
    )
    return df.with_columns(
        pl.exclude(config.GENE_ID_COLNAME).map_batches(
            lambda x: quantile_transform(x.to_frame(), **kwargs).flatten(),
            return_dtype=pl.Float64,
        )
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    set_max_resources(1, "8 GB", limit_polars=True)

    logger.info(f"Parsing {args.count_file.name}")
    count_df = parse_count_table(args.count_file)

    logger.info(f"Quantile normalising {args.count_file.name}")
    quantile_normalized_counts = quantile_normalise(count_df, args.target_distribution)

    export_parquet(quantile_normalized_counts, args.count_file, OUTFILE_SUFFIX)


if __name__ == "__main__":
    main()
