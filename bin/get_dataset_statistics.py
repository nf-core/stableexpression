#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

QUANT_NORM_SUFFIX = ".quant_norm.parquet"
DATASET_STATISTICS_SUFFIX = ".dataset_stats.csv"


ALLOWED_TARGET_DISTRIBUTIONS = ["normal", "uniform"]


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
    parser.add_argument(
        "--output", type=str, dest="outfile_name", required=True, help="Output file"
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


def compute_kolmogorov_smirnov_test_to_target_distribution(
    count_df: pd.DataFrame, target_distribution: str
) -> pd.Series:
    """Compute Kolmogorov-Smirnov test to target distribution."""

    if target_distribution == "normal":
        cum_distrib_function = stats.norm.cdf
    elif target_distribution == "uniform":
        cum_distrib_function = stats.uniform.cdf
    else:
        raise ValueError(f"Unknown target distribution: {target_distribution}")

    ks_tests = pd.Series(index=count_df.columns)
    for col in count_df.columns:
        ks = stats.ks_1samp(count_df[col], cum_distrib_function, nan_policy="omit")
        ks_tests[col] = ks.pvalue

    return ks_tests


def compute_dataset_statistics(
    count_df: pd.DataFrame, target_distribution: str
) -> pd.DataFrame:
    dataset_stats_df = count_df.describe()
    dataset_stats_df.loc["skewness"] = count_df.skew()
    # for each sample, test distance to target distribution
    ks_tests = compute_kolmogorov_smirnov_test_to_target_distribution(
        count_df, target_distribution
    )
    dataset_stats_df.loc[config.KS_TEST_COLNAME] = ks_tests
    return dataset_stats_df.T


def export_count_data(dataset_stats_df: pd.DataFrame, outfile_name: str):
    """Export dataset statistics to CSV files."""
    logger.info(f"Exporting dataset statistics counts to: {outfile_name}")
    dataset_stats_df.index.name = config.SAMPLE_COLNAME
    dataset_stats_df.to_csv(outfile_name, index=True, header=True)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_file = args.count_file

    logger.info(f"Computing dataset statistics for {count_file.name}")
    count_df = pd.read_parquet(count_file)
    count_df.set_index(config.GENE_ID_COLNAME, inplace=True)

    dataset_stats_df = compute_dataset_statistics(count_df, args.target_distribution)

    export_count_data(dataset_stats_df, args.outfile_name)


if __name__ == "__main__":
    main()
