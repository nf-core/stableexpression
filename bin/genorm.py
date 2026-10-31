#!/usr/bin/env python3

# Written by Olivier Coen. Highly optimised by Mistral. Released under the MIT license.
"""
genorm.py - Compute M-measure for each gene (scalable NumPy + Numba version).

Polars is used for I/O only (parquet in, csv out).

Key optimisation: the std of the pairwise log-ratio of two genes never
requires materialising the G*(G-1)/2 gene pairs. With a = log2(gene i) and
b = log2(gene j) across samples, and â, b̂ their row-centred profiles:

    Var0(a - b) = (|â|² + |b̂|² - 2·â·b̂) / S

The pairwise dot products â·b̂ form a Gram matrix obtained with a single BLAS
matmul per gene block, so the full (n_genes x n_genes) matrix is never
materialised. This drops memory from O(G²·S) pair rows to one (G x S) matrix
plus (block x G) temporaries, and lets BLAS do the heavy lifting.
"""

import argparse
import logging
import sys
from pathlib import Path

import config
from common import export_parquet, get_count_columns
import numpy as np
import polars as pl
from numba import njit, prange
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# number of rows processed per Gram block (tunes peak memory, not the result)
BLOCK_SIZE = 1024
EXPORT_NB_DIGITS = 3


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PARSE ARGS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--file",
        type=Path,
        dest="count_file",
        required=True,
        help="File where each row contains counts for two genes",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="outfile",
        required=True,
        help="M-measure output file",
    )
    return parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NUMBA-ACCELERATED KERNEL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@njit(parallel=True, cache=True)
def _m_sums_for_block(var_block: np.ndarray) -> np.ndarray:
    """
    var_block: (block, n_genes) matrix of Var0(a_i - b_j) / n_samples
    Returns, for each row i, the sum over j of sqrt(var) with negatives
    (numerical cancellation) treated as zero.
    """
    n_block, n_genes = var_block.shape
    out = np.empty(n_block, dtype=np.float64)
    for i in prange(n_block):
        s = 0.0
        for j in range(n_genes):
            v = var_block[i, j]
            if v > 0.0:
                s += np.sqrt(v)
        out[i] = s
    return out


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMPUTATION OF M-MEASURES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def compute_m_measures(L: np.ndarray, block_size: int = BLOCK_SIZE) -> np.ndarray:
    """
    L: (n_genes, n_samples) matrix of log2 counts.
    Returns the M-measure of each gene: the mean over all other genes of the
    std (ddof=0, across samples) of the pairwise log2 expression ratios.
    """
    n_genes, n_samples = L.shape

    if n_genes < 2:
        logger.error("Too few genes to compute M-measures")
        sys.exit(100)

    # row-centred profiles: std(log2(x_i) - log2(x_j)) only depends on these
    C = L - L.mean(axis=1, keepdims=True)  # (n_genes, n_samples)
    q = np.einsum("ij,ij->i", C, C)  # |â|^2 per gene, (n_genes,)

    m_measures = np.empty(n_genes, dtype=np.float64)

    for a in tqdm(range(0, n_genes, block_size), desc="geNorm blocks"):
        b = min(a + block_size, n_genes)

        # Gram block: dot products of gene rows a..b with all gene rows (BLAS)
        gram = C[a:b] @ C.T  # (block, n_genes)

        # Var0(log ratio) for every pair (i in block, j all genes)
        var = (q[a:b, None] + q[None, :] - 2.0 * gram) / n_samples

        # a gene paired with itself has zero variance (also guards the diagonal)
        var[np.arange(b - a), np.arange(a, b)] = 0.0

        # sum of stds over all partners for each gene of this block
        m_measures[a:b] = _m_sums_for_block(var)

    # each gene is paired with (n_genes - 1) partners
    m_measures /= n_genes - 1
    return m_measures


#####################################################
# MAIN
#####################################################

def main():
    args = parse_args()

    logger.info(f"Computing M-measure for {args.count_file}")

    count_lf = pl.scan_parquet(args.count_file)
    count_cols = get_count_columns(count_lf)

    # single collect of the count matrix
    count_df = count_lf.collect()
    genes = count_df[config.GENE_ID_COLNAME].to_list()
    X = (
        count_df
        .select(count_cols)
        .cast(pl.Float32)
        .fill_null(np.nan)
        .to_numpy()
    )  # (n_genes, n_samples), row order = genes

    if not np.isfinite(X).all() or (X < 0).any():
        logger.error(
            "Count matrix contains missing, non-finite or negative values; "
            "geNorm requires complete positive counts (impute first)"
        )
        sys.exit(100)

    # replace 0 values by a very small count value
    X[X == 0] = 1e-8
    L = np.log2(X)

    m_measures = compute_m_measures(L)

    m_df = (
        pl.DataFrame(
            {
                config.GENE_ID_COLNAME: genes,
                config.GENORM_M_MEASURE_COLNAME: m_measures,
            }
        )
        .sort(config.GENE_ID_COLNAME)
    )

    # limit nb of digits to 5 to avoid discrepancies between docker / apptainer and conda runs
    export_parquet(m_df, args.outfile, EXPORT_NB_DIGITS)

    logger.info("Done")


if __name__ == "__main__":
    main()
