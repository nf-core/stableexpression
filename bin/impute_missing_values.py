#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import export_parquet, parse_count_table
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer, KNNImputer

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


def separate_genes_with_high_number_of_zeros(
    df: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Separate genes with high number of zeros from the rest of the genes.
    """
    max_number_of_zeros = THRESHOLD_RATIO_ZEROS * len(get_count_columns(df))
    to_keep = (
        pl.sum_horizontal(pl.exclude(config.GENE_ID_COLNAME).eq(0))
        <= max_number_of_zeros
    )
    return df.filter(to_keep), df.filter(~to_keep)


def replace_nulls_with_row_mean(df: pl.DataFrame):
    return (
        df.with_columns(
            pl.mean_horizontal(
                pl.exclude(config.GENE_ID_COLNAME), ignore_nulls=True
            ).alias("row_mean")
        )
        .with_columns(
            [pl.col(col).fill_null(pl.col("row_mean")) for col in get_count_columns(df)]
        )
        .drop("row_mean")
    )


def apply_knn_imputer(df: pl.DataFrame):
    logger.info("Applying KNN imputation")
    imputer = KNNImputer(n_neighbors=N_NEIGHBORS, weights="distance")
    # Convert to numpy, impute, then convert back
    count_matrix = df.select(get_count_columns(df)).to_numpy()
    imputed_array = imputer.fit_transform(count_matrix)
    return df.with_columns(pl.DataFrame(imputed_array, schema=get_count_columns(df)))


def apply_iterative_imputer(df: pl.DataFrame):
    logger.info("Applying iterative imputation")
    imputer = IterativeImputer(
        max_iter=MAX_ITERATIONS,
        n_nearest_features=N_NEAREST_FEATURES,
        random_state=0,
        initial_strategy="mean",
    )
    # Convert to numpy, impute, then convert back
    count_matrix = df.select(get_count_columns(df)).to_numpy()
    imputed_array = imputer.fit_transform(count_matrix)
    return df.with_columns(pl.DataFrame(imputed_array, schema=get_count_columns(df)))


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_file = args.count_file

    logger.info(f"Parsing {count_file.name}")
    df = parse_count_table(count_file)

    # logger.info("Separating genes with high number of zeros")
    # df, high_zero_genes_df = separate_genes_with_high_number_of_zeros(count_df)

    if args.imputer == "iterative":
        df = apply_iterative_imputer(df)
    elif args.imputer == "knn":
        df = apply_knn_imputer(df)
    elif args.imputer == "gene_mean":
        df = replace_nulls_with_row_mean(df)

    export_parquet(df, count_file, OUTFILE_SUFFIX)


if __name__ == "__main__":
    main()
