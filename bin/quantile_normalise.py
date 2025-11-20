#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import pandas as pd
from sklearn.preprocessing import QuantileTransformer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

QUANT_NORM_SUFFIX = ".quant_norm.parquet"

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
    return parser.parse_args()


def quantile_normalise(data: pd.DataFrame, target_distribution: str):
    """
    Quantile normalize a data matrix based on a target distribution.
    """
    transformer = QuantileTransformer(
        n_quantiles=N_QUANTILES, output_distribution=target_distribution, subsample=None
    )

    normalised_data = pd.DataFrame(index=data.index, columns=data.columns)
    for col in data.columns:
        normalised_data[col] = transformer.fit_transform(data[col].to_frame())

    return normalised_data


def export_count_data(quantile_normalized_counts: pd.DataFrame, count_file: Path):
    """Export gene expression data to CSV files."""
    outfilename = count_file.name.replace(".csv", QUANT_NORM_SUFFIX)
    logger.info(f"Exporting quantile normalised counts to: {outfilename}")
    quantile_normalized_counts.reset_index().to_parquet(outfilename)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_file = args.count_file

    logger.info(f"Quantile normalising {count_file.name}")
    # count_df = pd.read_parquet(count_file)
    # count_df.set_index(config.GENE_ID_COLNAME, inplace=True)
    count_df = pd.read_csv(count_file, index_col=0)
    count_df.index.name = config.GENE_ID_COLNAME

    quantile_normalized_counts = quantile_normalise(count_df, args.target_distribution)

    export_count_data(quantile_normalized_counts, count_file)


if __name__ == "__main__":
    main()
