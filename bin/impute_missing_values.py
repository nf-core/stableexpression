#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
This script imputes missing values in a count table using KNN imputation.
It takes as inputs:
    - a dataframe whose columns are samples and rows are genes.
    - an output file.
Missing value imputation is performed using fknni (FAISS-backed KNN imputer).
To scale to large datasets, the genes (features) are split into random blocks
that are imputed independently. All samples are kept in each block, so the
pool of KNN neighbours (samples) is never restricted.
"""

import argparse
import logging
from pathlib import Path

import config
import numpy as np
import polars as pl
from common import export_parquet, get_count_columns
from fknni import FastKNNImputer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

N_NEIGHBOURS = 5

GENE_BLOCK_SIZE = 200 # sot that 20000 genes are split into 100 blocks

RANDOM_STATE = 42


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
    Rows are samples and columns are genes so that the imputer uses gene
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


def make_gene_blocks(n_genes: int, block_size: int, rng: np.random.Generator) -> list[np.ndarray]:
    """
    Split gene indices into random blocks of at most `block_size` genes.

    The assignment is random (shuffled) so that each block is representative of
    the full gene set: contiguous slices could group correlated genes and bias
    the per-block distances.

    Returns the list of gene index arrays, in the original gene order.
    """
    shuffled = rng.permutation(n_genes)
    blocks = [shuffled[i : i + block_size] for i in range(0, n_genes, block_size)]
    # keep original gene order within each block so reassembly is trivial
    return [np.sort(block) for block in blocks]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
    args = parse_args()

    logger.info(f"Parsing {args.count_file.name}")
    lf = pl.scan_parquet(args.count_file)

    logger.info("Applying FAISS KNN imputation")
    count_cols, count_matrix = get_count_matrix(lf)
    n_samples, n_genes = count_matrix.shape

    rng = np.random.default_rng(RANDOM_STATE)
    gene_blocks = make_gene_blocks(n_genes, GENE_BLOCK_SIZE, rng)
    logger.info(
        f"Imputing {n_genes} genes in {len(gene_blocks)} blocks of " +
        f"~{GENE_BLOCK_SIZE} genes ({n_samples} samples per block)"
    )

    imputer = FastKNNImputer(n_neighbors=N_NEIGHBOURS)

    # impute each gene block independently, in place
    for i, gene_idx in enumerate(gene_blocks):
        logger.info(f"Imputing block {i + 1}/{len(gene_blocks)} ({len(gene_idx)} genes)")
        block = count_matrix[:, gene_idx]
        count_matrix[:, gene_idx] = imputer.fit_transform(block)

    df = pl.concat(
        [
            lf.select(config.GENE_ID_COLNAME).collect(),
            pl.DataFrame(count_matrix.T, schema=count_cols),  # shape (n_genes, n_samples)
        ],
        how="horizontal",
    )

    export_parquet(df, args.outfile)

    logger.info("Done")


if __name__ == "__main__":
    main()
