#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import export_parquet, parse_count_table
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".imputed.parquet"

THRESHOLD_RATIO_ZEROS = 0.9

# KNN
N_NEIGHBORS = 10

# ITERATIVE
MAX_ITERATIONS = 10
N_NEAREST_FEATURES = 100

IMPUTERS = ["knn", "iterative", "gene_mean"]


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Perform KNN imputation on count data")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument("--imputer", choices=IMPUTERS, required=True, dest="imputer")
    return parser.parse_args()


def get_count_columns(df: pl.DataFrame):
    return df.select(pl.exclude(config.GENE_ID_COLNAME)).columns


def apply_imputer(df: pl.DataFrame, imputer):
    # convert to numpy, impute, then convert back
    count_matrix = df.select(get_count_columns(df)).to_numpy()
    imputed_array = imputer.fit_transform(count_matrix)
    return df.with_columns(pl.DataFrame(imputed_array, schema=get_count_columns(df)))


def apply_simle_imputer(df: pl.DataFrame):
    imputer = SimpleImputer()
    return apply_imputer(df, imputer)


def apply_knn_imputer(df: pl.DataFrame) -> pl.DataFrame:
    imputer = KNNImputer(n_neighbors=N_NEIGHBORS, weights="distance")
    return apply_imputer(df, imputer)


def apply_iterative_imputer(df: pl.DataFrame) -> pl.DataFrame:
    imputer = IterativeImputer(
        max_iter=MAX_ITERATIONS,
        n_nearest_features=N_NEAREST_FEATURES,
        random_state=0,
        initial_strategy="mean",
        min_value=0,
        max_value=1,
        imputation_order="random",
        verbose=1,
    )
    return apply_imputer(df, imputer)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info(f"Parsing {args.count_file.name}")
    df = parse_count_table(args.count_file)

    # logger.info("Separating genes with high number of zeros")
    # df, high_zero_genes_df = separate_genes_with_high_number_of_zeros(count_df)

    if args.imputer == "iterative":
        logger.info("Applying iterative imputation")
        df = apply_iterative_imputer(df)
    elif args.imputer == "knn":
        logger.info("Applying KNN imputation")
        df = apply_knn_imputer(df)
    elif args.imputer == "gene_mean":
        logger.info("Applying simple imputation")
        df = apply_simle_imputer(df)

    export_parquet(df, args.count_file, OUTFILE_SUFFIX)

    logger.info("Done")


if __name__ == "__main__":
    main()
