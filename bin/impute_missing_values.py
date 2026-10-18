#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
This script imputes missing values in a count table using a specified imputation method.
It takes as inputs:
    - a dataframe whose columns are samples and rows are genes.
    - the name of the imputation method to use.
Missing value imputation is performed using scikit-learn.
Prior to imputation, a Mini Batch K-means algorithm is used to cluster the genes.
"""

import argparse
import logging
from pathlib import Path

import config
import numpy as np
import polars as pl
from common import export_parquet, get_count_columns, get_nb_rows
from fknni import FastKNNImputer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

MIN_NEIGHBOURS = 10
MAX_NEIGHBOURS = 30

FACTOR_N_SAMPLES_TO_K = 2


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CLI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_args():
    parser = argparse.ArgumentParser(description="Perform KNN imputation on count data")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--out", type=Path, dest="outfile", required=True, help="Output file"
    )
    return parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MATRIX HELPERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_count_matrix(lf: pl.LazyFrame) -> tuple[list[str], np.ndarray]:
    """
    Extract the count matrix as float32 numpy array of shape (n_samples, n_genes).
    Rows are samples and columns are genes so that the imputers use gene
    expression values as features. The transposition is performed on the numpy
    array because it is zero-copy (it only changes memory strides).

    float32 halves the memory footprint compared to float64, which is required
    for large datasets (e.g. 25k genes x 35k samples ~ 3.5 GB in float32).
    """
    count_cols = get_count_columns(lf)
    count_matrix = (
        lf.select(count_cols)
        .collect()
        .cast(pl.Float32)
        .to_numpy()
        .T
    )  # shape (n_samples, n_genes)
    return count_cols, count_matrix


def matrix_to_dataframe(
    lf: pl.LazyFrame, count_cols: list[str], imputed_array: np.ndarray
) -> pl.DataFrame:
    """Convert the imputed array (n_samples, n_genes) back to a dataframe (genes x samples)."""
    return pl.concat(
        [
            lf.select(config.GENE_ID_COLNAME).collect(),
            pl.DataFrame(imputed_array.T, schema=count_cols),  # shape (n_genes, n_samples)
        ],
        how="horizontal",
    )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPUTERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_number_of_neighbours(lf: pl.LazyFrame, k_min: int, k_max: int) -> int:
    """
    Returns the number of neighbours to use for KNN-imputation based on the
    number of samples and missing values. The number of neighbours is
    determined by the square root of the number of samples, scaled by a
    factor to account for missing values. It is clipped to [k_min, k_max].
    """
    # proportion of missing values in the dataframe
    nb_missing_values = (
        lf.select(pl.exclude(config.GENE_ID_COLNAME))
        .null_count()
        .collect()
        .sum_horizontal()
        .item()
    )
    n_genes = get_nb_rows(lf)
    n_samples = len(get_count_columns(lf))
    missing_fraction = nb_missing_values / (n_genes * n_samples)
    # return k as a function of the number of samples and missing values, with bounds
    return int(
        np.clip(
            int(np.sqrt(n_samples) / FACTOR_N_SAMPLES_TO_K * (1 + missing_fraction)),
            k_min,
            k_max,
        )
    )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
    args = parse_args()

    logger.info(f"Parsing {args.count_file.name}")
    lf = pl.scan_parquet(args.count_file)

    logger.info("Applying FAISS KNN imputation")
    count_cols, count_matrix = get_count_matrix(lf)

    n_neighbours = get_number_of_neighbours(lf, MIN_NEIGHBOURS, MAX_NEIGHBOURS)
    logger.info(f"Using {n_neighbours} neighbours")

    # KNN imputation on the full matrix using fknni (FAISS-backed, multithreaded)
    imputer = FastKNNImputer(n_neighbors=n_neighbours)
    imputed_array = imputer.fit_transform(count_matrix)  # shape (n_samples, n_genes)

    df = matrix_to_dataframe(lf, count_cols, imputed_array)

    export_parquet(df, args.outfile)

    logger.info("Done")


if __name__ == "__main__":
    main()
