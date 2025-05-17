#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
from pathlib import Path
from scipy import stats
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

QUANT_NORM_SUFFIX = ".quant_norm.parquet"
DATASET_STATISTICS_SUFFIX = ".dataset_stats.csv"

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
SAMPLE_COLNAME = "sample"
KS_TEST_COLNAME = "kolmogorov_smirnov_to_uniform_dist_pvalue"


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
    return parser.parse_args()


def compute_kolmogorov_smirnov_test_to_uniform_distribution(count_df: pd.DataFrame):
    """Compute Kolmogorov-Smirnov test to uniform distribution."""
    ks_tests = pd.Series(index=count_df.columns)
    for col in count_df.columns:
        ks = stats.ks_1samp(count_df[col], stats.uniform.cdf, nan_policy="omit")
        ks_tests[col] = ks.pvalue
    return ks_tests


def compute_dataset_statistics(count_df: pd.DataFrame):
    dataset_stats_df = count_df.describe()
    dataset_stats_df.loc["skewness"] = count_df.skew()
    # for each sample, test distance to uniform distribution
    ks_tests = compute_kolmogorov_smirnov_test_to_uniform_distribution(count_df)
    dataset_stats_df.loc[KS_TEST_COLNAME] = ks_tests
    return dataset_stats_df.T


def export_count_data(dataset_stats_df: pd.DataFrame, outfile_name: str):
    """Export dataset statistics to CSV files."""
    logger.info(f"Exporting dataset statistics counts to: {outfile_name}")
    dataset_stats_df.index.name = SAMPLE_COLNAME
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
    count_df.set_index(ENSEMBL_GENE_ID_COLNAME, inplace=True)

    dataset_stats_df = compute_dataset_statistics(count_df)

    export_count_data(dataset_stats_df, args.outfile_name)


if __name__ == "__main__":
    main()
