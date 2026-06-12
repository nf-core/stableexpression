#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import numpy as np
import polars as pl
from common import export_parquet, parse_count_table
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.cluster import MiniBatchKMeans

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".imputed.parquet"

IMPUTERS = ["knn", "iterative", "gene_mean"]

MAX_ITER = 10
N_NEAREST_FEATURES = 50

MIN_NEIGHBOURS = 10
MAX_NEIGHBOURS = 50

FACTOR_N_SAMPLES_TO_K = 2

BASE_N_CLUSTERS = 10

KEEP_EMPTY_FEATURES = True


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Perform KNN imputation on count data")
    parser.add_argument(
        "--counts", 
        type=Path, 
        dest="count_file", 
        required=True, 
        help="Count file"
    )
    parser.add_argument(
        "--imputer", 
        choices=IMPUTERS, 
        required=True, 
        dest="imputer",
        help="Imputer to use"
    )    
    return parser.parse_args()


def get_count_columns(df: pl.DataFrame):
    return df.select(pl.exclude(config.GENE_ID_COLNAME)).columns


def apply_imputer(df: pl.DataFrame, imputer):
    # convert to numpy, impute, then convert back
    # we transpose the matrix so that each row is a sample and each column is a gene
    # this is so that the imputer can use the gene expression values as features
    # we transpose once it is a numpy array because numpty transposition is zero-copy transpose (just changes memory strides)
    # the to_numpy method should be zero copy too, since all count columns have a float dtype
    count_matrix = df.select(get_count_columns(df)).to_numpy().T # (n_samples, n_genes)
    imputed_array = imputer.fit_transform(count_matrix)   
    return df.with_columns(pl.DataFrame(imputed_array.T, schema=get_count_columns(df)))


def apply_simle_imputer(df: pl.DataFrame):
    imputer = SimpleImputer(
        strategy="median",
        copy=False,
        keep_empty_features=KEEP_EMPTY_FEATURES,
    )
    return apply_imputer(df, imputer)


def get_number_of_neighbours(df: pl.DataFrame, k_min: int, k_max: int) -> int:
    """
    Returns the number of neighbours to use for KNN-imputation based on the number of samples and missing values.
    
    Parameters
    ----------
    df : pl.DataFrame
        The input dataframe.
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
        df.select(pl.exclude(config.GENE_ID_COLNAME))
        .null_count()
        .sum_horizontal()
        .item()
    )
    n_genes = len(df)
    n_samples = df.select(pl.exclude(config.GENE_ID_COLNAME)).width
    missing_fraction = nb_missing_values / (n_genes * n_samples)
    # return k as a function of the number of samples and missing values, with bounds
    return int(np.clip(
        int(np.sqrt(n_samples) / FACTOR_N_SAMPLES_TO_K * (1 + missing_fraction)),
        k_min, k_max
    ))


def cluster_dataframe(df: pl.DataFrame, n_clusters: int) -> pl.Series:
    """
    Cluster the dataframe using MiniBatchKMeans.
    
    Returns
    -------
    pl.Series
        Cluster labels for each gene.
    """
    # replace missing values with mean values accross samples (for clustering only)
    df_for_clustering = (
        df.select(pl.exclude(config.GENE_ID_COLNAME))
        .with_columns(
            mean_expression=pl.mean_horizontal(pl.all())
        )
        .fill_null(pl.col("mean_expression"))
        .drop("mean_expression")
    )
    
    # cluster — MiniBatchKMeans scales well
    kmeans = MiniBatchKMeans(
        n_clusters=n_clusters,
        random_state=42,
    )
    # return cluster labels (gene per gene)
    return pl.Series(kmeans.fit_predict(df_for_clustering))


def get_number_of_clusters(df: pl.DataFrame) -> int:
    """
    Returns the number of clusters to use for KNN-imputation based on the number of samples and missing values.
    In the very rare case where the number of samples is less than the base number of clusters, the number of clusters is set to the number of samples.
    This is especially useful for some specific test cases where the number of samples is very small.
    Parameters
    ----------
    df : pl.DataFrame
        The input dataframe.
    
    Returns
    -------
    int
        The number of clusters to use for KNN-imputation.
    """
    return min(BASE_N_CLUSTERS, df.height)


def apply_knn_imputer(df: pl.DataFrame) -> pl.DataFrame:
    n_neighbours = get_number_of_neighbours(df, MIN_NEIGHBOURS, MAX_NEIGHBOURS)
    logger.info(f"Using {n_neighbours} neighbours")
    imputer = KNNImputer(
        n_neighbors=n_neighbours, 
        weights="distance",
        copy=False,
        keep_empty_features=KEEP_EMPTY_FEATURES,
    )

    n_clusters = get_number_of_clusters(df)
    logger.info(f"Making {n_clusters} clusters using MiniBatchKMeans")
    labels = cluster_dataframe(df, n_clusters=n_clusters)

    unique_labels = labels.unique()
    for label in unique_labels:
        cluster_mask = labels == label
        cluster_df = df.filter(cluster_mask)
        logger.info(f"Imputing cluster {label} ({cluster_df.shape[0]} genes)")
        imputed_cluster_df = apply_imputer(cluster_df, imputer)
        imputed_cluster_df.write_parquet(f"imputed_cluster_{label}.parquet")

    del df
    return pl.concat([pl.read_parquet(f"imputed_cluster_{label}.parquet") for label in unique_labels])


def apply_iterative_imputer(df: pl.DataFrame) -> pl.DataFrame:
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

    if args.imputer == "knn":
        logger.info("Applying KNN imputation")
        df = apply_knn_imputer(df)
    elif args.imputer == "iterative":
        logger.info("Applying iterative imputation")
        df = apply_iterative_imputer(df)
    elif args.imputer == "gene_mean":
        logger.info("Applying simple imputation")
        df = apply_simle_imputer(df)

    export_parquet(df, args.count_file, OUTFILE_SUFFIX)

    logger.info("Done")


if __name__ == "__main__":
    main()
