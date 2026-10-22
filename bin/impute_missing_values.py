#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
This script imputes missing values in a count table using KNN imputation.
It takes as inputs:
    - a dataframe whose columns are samples and rows are genes.
    - an output file.
Missing value imputation is performed using scikit-learn's KNNImputer.
To scale to large datasets, the samples are split into random groups of at
most MAX_GROUP_SIZE samples, and each group is imputed independently.
Multiple passes are performed with different random splits: imputations from
all passes are averaged per cell, which reduces the variance introduced by
restricting the neighbour pool to a single group. Groups are imputed in
parallel with joblib.
"""

import argparse
import logging
from pathlib import Path
from math import sqrt

import config
import numpy as np
import polars as pl
from common import export_parquet, get_count_columns
from joblib import Parallel, delayed
from sklearn.impute import KNNImputer

PROCESS_LOGGING_FORMAT = "%(levelname)s:%(processName)s: %(message)s"
N_NEIGHBOURS = 5
MAX_GROUP_SIZE = 1000
MAX_N_PASSES = 5
RANDOM_STATE = 42

logging.basicConfig(level=logging.INFO, format=PROCESS_LOGGING_FORMAT, force=True)
logger = logging.getLogger(__name__)


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
    parser.add_argument(
        "--threads",
        type=int,
        required=True,
        help="Number of parallel jobs for joblib",
    )
    return parser.parse_args()


def has_missing_values(lf: pl.LazyFrame) -> bool:
    """
    Check if there are any missing values in the input counts.
    """
    return lf.select(
        pl.any_horizontal(pl.exclude(config.GENE_ID_COLNAME).is_null()).any()
    ).collect().item()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MATRIX HELPERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_count_matrix(lf: pl.LazyFrame) -> tuple[list[str], np.ndarray]:
    """
    Extract the count matrix as float32 numpy array of shape (n_samples, n_genes).
    Rows are samples and columns are genes so that the imputer uses gene
    expression values as features. The transposition is performed on the numpy
    array because it is zero-copy (it only changes memory strides).

    float32 halves the memory footprint compared to float64, which is required
    for large datasets (e.g. 25k genes x 35k samples ~ 3.5 GB in float32).
    """
    count_cols = get_count_columns(lf)
    count_matrix = (
        lf.select(count_cols)
        .fill_null(np.nan) # convert all null values to NaN
        .collect()
        .cast(pl.Float32)
        .to_numpy()  # replaces null values by NaNs
        .T
    )  # shape (n_samples, n_genes)
    return count_cols, count_matrix


def get_global_gene_means(count_matrix: np.ndarray) -> np.ndarray:
    """
    Per-gene mean computed over all samples, ignoring NaNs.
    Genes with no observed value at all are mapped to 0.0 (safety net).
    """
    with np.errstate(invalid="ignore"):
        gene_means = np.nanmean(count_matrix, axis=0)
    return np.where(np.isnan(gene_means), 0.0, gene_means)


def make_random_groups(
    n_samples: int, group_size: int, rng: np.random.Generator
) -> list[np.ndarray]:
    """
    Randomly split sample indices into groups of at most `group_size` samples.
    A new call with a different rng state yields a different split.
    """
    shuffled = rng.permutation(n_samples)
    return [shuffled[i : i + group_size] for i in range(0, n_samples, group_size)]


def setup_worker_logging():
    """
    Set up logging for worker processes.
    """
    logging.basicConfig(level=logging.INFO, format=PROCESS_LOGGING_FORMAT, force=True)


def get_nb_passes(n_samples: int):
    """
    Compute the number of passes over the data.
    The number of passes is the square root of the number of samples divided by
    the maximum group size.
    It is clipped between 1 and MAX_N_PASSES.
    """
    n_passes = sqrt(n_samples // MAX_GROUP_SIZE)
    return int(np.clip(n_passes, 1 , MAX_N_PASSES))


def process_single_group(
    sample_idx: np.ndarray,
    sub_matrix: np.ndarray,
    group_nan_mask: np.ndarray,
    group_id: int,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray] | None:
    """
    Worker function executed in parallel processes (loky backend).
    Processes a single chunk of samples.

    Returns a tuple (sample_idx, valid_genes_idx, filled, sub_nan_mask)
    where filled contains the imputed values (0.0 where the cell was observed
    or not imputed), or None if the chunk is uncomputable (all genes fully
    missing in this chunk — other passes usually cover it).
    """

    fully_missing = group_nan_mask.all(axis=0)
    if fully_missing.all():
        logger.warning(f"Group {group_id}: fully missing")
        # all genes have only missing values in this chunk: nothing to impute
        return None

    imputer = KNNImputer(
        n_neighbors=N_NEIGHBOURS,
        weights="distance",
        keep_empty_features=True,
    )

    if fully_missing.any():
        ratio_missing = group_nan_mask.sum() / group_nan_mask.size
        logger.info(f"Group {group_id}: {ratio_missing:.1%} genes with missing values")
        # some genes have only missing values: impute on the valid genes only
        valid_genes_idx = np.where(~fully_missing)[0]
        sub_matrix_input = sub_matrix[:, valid_genes_idx]
        sub_nan_mask = group_nan_mask[:, valid_genes_idx]

        sub_imputed = imputer.fit_transform(sub_matrix_input)
        filled = np.where(sub_nan_mask, sub_imputed, 0.0)

        return sample_idx, valid_genes_idx, filled, sub_nan_mask

    # all genes have at least one non-missing value
    else:
        logger.info(f"Group {group_id}: no missing value")
        sub_imputed = imputer.fit_transform(sub_matrix)
        filled = np.where(group_nan_mask, sub_imputed, 0.0)

        return sample_idx, None, filled, group_nan_mask


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MULTI-PASS IMPUTATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def multi_pass_impute(count_matrix: np.ndarray, n_passes: int, nthreads: int, rng: np.random.Generator) -> np.ndarray:
    """
    Impute the matrix with possibly multiple passes over random sample groups.
    On most dataframes, 1 pass is enough to cover all missing values.
    Multiple passes are used when count matrix with fully missing values are encountered.
    Parallelizes the chunks over the available CPUs with joblib (loky).
    """
    nan_mask = np.isnan(count_matrix)

    # store imputed values (float64 to avoid precision loss)
    imputed_sum = np.zeros(count_matrix.shape, dtype=np.float64)
    imputed_count = np.zeros(count_matrix.shape, dtype=np.int32)

    for pass_nb in range(1, n_passes + 1):
        groups = make_random_groups(count_matrix.shape[0], MAX_GROUP_SIZE, rng)
        logger.info(f"Pass {pass_nb}/{n_passes}: {len(groups)} groups")

        results = Parallel(
            n_jobs=nthreads,
            backend="loky",
            initializer=setup_worker_logging,
        )(
            delayed(process_single_group)(
                sample_idx,
                count_matrix[sample_idx],
                nan_mask[sample_idx],
                group_id,
            )
            for group_id, sample_idx in enumerate(groups)
        )

        # safe sequential aggregation of the results from the different workers
        for res in results:
            if res is None:
                continue

            sample_idx, valid_genes_idx, filled, sub_mask = res

            if valid_genes_idx is not None:
                imputed_sum[sample_idx[:, None], valid_genes_idx] += filled
                imputed_count[sample_idx[:, None], valid_genes_idx] += sub_mask
            else:
                imputed_sum[sample_idx] += filled
                imputed_count[sample_idx] += sub_mask

        n_remaining = int((nan_mask & (imputed_count == 0)).sum())
        logger.info(
            f"Pass {pass_nb}/{n_passes} done: "
            f"{100 * (imputed_count > 0).sum() / max(nan_mask.sum(), 1):.1%} of "
            f"missing cells imputed at least once ({n_remaining} remaining)"
        )

    # means of imputations (initial values stay unchanged)
    result = count_matrix.astype(np.float64, copy=True)
    has_updates = nan_mask & (imputed_count > 0)
    result[has_updates] = imputed_sum[has_updates] / imputed_count[has_updates]

    return result.astype(np.float32)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
    args = parse_args()

    logger.info(f"Parsing {args.count_file.name}")
    lf = pl.scan_parquet(args.count_file)

    if has_missing_values(lf):

        logger.info("Applying KNN imputation")
        count_cols, count_matrix = get_count_matrix(lf)
        n_samples, n_genes = count_matrix.shape
        logger.info(f"Found {n_samples} samples and {n_genes} genes")

        # pristine gene means, computed before any imputation
        gene_means = get_global_gene_means(count_matrix)

        n_passes = get_nb_passes(n_samples)
        logger.info(f"Will perform {n_passes} passes")

        rng = np.random.default_rng(RANDOM_STATE)
        count_matrix = multi_pass_impute(count_matrix, n_passes, args.threads, rng)

        # if some values are still missing after all passes
        still_missing = np.isnan(count_matrix)
        if still_missing.any():
            logger.info(
                f"{still_missing.sum()} cells still missing after all passes: "
                "filling with pristine global gene means"
            )
            count_matrix[still_missing] = np.take(gene_means, np.where(still_missing)[1])

        df = pl.concat(
            [
                lf.select(config.GENE_ID_COLNAME).collect(),
                pl.DataFrame(count_matrix.T, schema=count_cols),  # shape (n_genes, n_samples)
            ],
            how="horizontal",
        )

        logger.info(f"Exporting imputed counts to {args.outfile}")
        export_parquet(df, args.outfile)

    else:
        logger.info(f"No missing values found. Exporting counts to {args.outfile}")
        lf.sink_parquet(args.outfile)

    logger.info("Done")


if __name__ == "__main__":
    main()
