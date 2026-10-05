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
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.cluster import MiniBatchKMeans

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

IMPUTERS = ["knn", "iterative", "gene_mean"]

MAX_ITER = 10
N_NEAREST_FEATURES = 50

MIN_NEIGHBOURS = 10
MAX_NEIGHBOURS = 50

FACTOR_N_SAMPLES_TO_K = 2

BASE_N_CLUSTERS = 10

KEEP_EMPTY_FEATURES = True


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CLI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_args():
    parser = argparse.ArgumentParser(description="Perform KNN imputation on count data")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--imputer",
        choices=IMPUTERS,
        required=True,
        dest="imputer",
        help="Imputer to use",
    )
    parser.add_argument(
        "--out", type=Path, dest="outfile", required=True, help="Output file"
    )
    return parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CLUSTERING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_number_of_neighbours(lf: pl.LazyFrame, k_min: int, k_max: int) -> int:
    """
    Returns the number of neighbours to use for KNN-imputation based on the number of samples and missing values.
    The number of neighbours is determined by the square root of the number of samples, scaled by a factor to account for missing values.
    It is clipped to the range [k_min, k_max] to ensure the number of neighbours is within the specified bounds.

    Parameters
    ----------
    lf: pl.LazyFrame
        The input lazyframe.
    k_min : int
        The minimum number of neighbours.
    k_max : int
        The maximum number of neighbours.

    Returns
    -------
    int
        The number of neighbours to use for KNN-imputation.
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


def cluster_dataframe(lf: pl.LazyFrame, n_clusters: int) -> pl.Series:
    """
    Cluster the dataframe by similarity between genes using MiniBatchKMeans.

    Parameters
    ----------
    lf : pl.LazyFrame
        The input lazyframe.
    n_clusters : int
        The number of clusters to use.

    Returns
    -------
    pl.Series
        Cluster labels for each gene.
    """
    # replace missing values with mean values accross samples (for clustering only)
    df_for_clustering = (
        lf.select(pl.exclude(config.GENE_ID_COLNAME))
        .with_columns(mean_expression=pl.mean_horizontal(pl.all()))
        .fill_null(pl.col("mean_expression"))
        .drop("mean_expression")
        .collect()
    )

    # cluster — MiniBatchKMeans scales well
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=42,
    )
    # return cluster labels (gene per gene)
    return pl.Series(kmeans.fit_predict(df_for_clustering))


def get_number_of_clusters(lf: pl.LazyFrame) -> int:
    """
    Returns the number of clusters to use for KNN-imputation based on the number of samples and missing values.
    In the very rare case where the number of samples is less than the base number of clusters, the number of clusters is set to the number of samples.
    This is especially useful for some specific test cases where the number of samples is very small.
    Parameters
    ----------
    lf: pl.LazyFrame
        The input lazyframe.

    Returns
    -------
    int
        The number of clusters to use for KNN-imputation.
    """
    return min(BASE_N_CLUSTERS, get_nb_rows(lf))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPUTERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def apply_imputer(lf: pl.LazyFrame, imputer):
    # convert to numpy, impute, then convert back
    # we transpose the matrix so that each row is a sample and each column is a gene
    # this is so that the imputer can use the gene expression values as features
    # we transpose once it is a numpy array because numpty transposition is zero-copy transpose (just changes memory strides)
    # the to_numpy method should be zero copy too, since all count columns have a float dtype
    count_cols = get_count_columns(lf)
    count_matrix = lf.select(count_cols).collect().to_numpy().T  # shape (n_samples, n_genes)
    imputed_array = imputer.fit_transform(count_matrix) # shape (n_samples, n_genes)
    return pl.concat(
        [
            lf.select(config.GENE_ID_COLNAME).collect(),
            pl.DataFrame(imputed_array.T, schema=count_cols) # shape (n_genes, n_samples)
        ],
        how='horizontal'
    )


def apply_knn_imputer(lf: pl.LazyFrame) -> pl.DataFrame:
    n_neighbours = get_number_of_neighbours(lf, MIN_NEIGHBOURS, MAX_NEIGHBOURS)
    logger.info(f"Using {n_neighbours} neighbours")
    imputer = KNNImputer(
        n_neighbors=n_neighbours,
        weights="distance",
        copy=False,
        keep_empty_features=KEEP_EMPTY_FEATURES,
    )

    n_clusters = get_number_of_clusters(lf)
    logger.info(f"Making {n_clusters} clusters using MiniBatchKMeans")
    labels = cluster_dataframe(lf, n_clusters=n_clusters)

    unique_labels = labels.unique()
    label_parquet_files = []
    for label in unique_labels:
        cluster_mask = labels == label
        cluster_lf = lf.filter(cluster_mask)
        logger.info(f"Imputing cluster {label} ({get_nb_rows(cluster_lf)} genes)")
        imputed_cluster_df = apply_imputer(cluster_lf, imputer)
        file = f"imputed_cluster_{label}.parquet"
        imputed_cluster_df.write_parquet(file)
        label_parquet_files.append(file)

    df = pl.concat([pl.read_parquet(file) for file in label_parquet_files])
    # removing intermediate parquet files
    for file in label_parquet_files:
        Path(file).unlink()
    return df


def apply_iterative_imputer(lf: pl.LazyFrame) -> pl.DataFrame:
    imputer = IterativeImputer(
        max_iter=MAX_ITER,
        sample_posterior=True,
        n_nearest_features=N_NEAREST_FEATURES,
        random_state=42,
        initial_strategy="mean",
        min_value=0,
        max_value=1,
        imputation_order="random",
        keep_empty_features=KEEP_EMPTY_FEATURES,
        verbose=1,
    )
    return apply_imputer(lf, imputer)


def apply_simle_imputer(lf: pl.LazyFrame):
    imputer = SimpleImputer(
        strategy="median",
        copy=False,
        keep_empty_features=KEEP_EMPTY_FEATURES,
    )
    return apply_imputer(lf, imputer)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
    args = parse_args()

    logger.info(f"Parsing {args.count_file.name}")
    lf = pl.scan_parquet(args.count_file)

    if args.imputer == "knn":
        logger.info("Applying KNN imputation")
        df = apply_knn_imputer(lf)
    elif args.imputer == "iterative":
        logger.info("Applying iterative imputation")
        df = apply_iterative_imputer(lf)
    elif args.imputer == "gene_mean":
        logger.info("Applying simple imputation")
        df = apply_simle_imputer(lf)
    else:
        raise ValueError(f"Unknown imputer: {args.imputer}")

    export_parquet(df, args.outfile)

    logger.info("Done")


if __name__ == "__main__":
    main()
