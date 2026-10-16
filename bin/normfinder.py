#!/usr/bin/env python3

# NormFinder stability values — scalable NumPy + Numba implementation.
# Polars is used for I/O only (parquet/csv in, csv out).
# All linear algebra runs on a single dense (n_genes x n_samples) float32 matrix,
# with NaN as the canonical "missing value" marker.
#
# Written originally by Olivier Coen. Highly optimised by Mistral. Released under the MIT license.

import argparse
import logging
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import config
import numpy as np
import polars as pl
from common import write_csv_with_floats
from numba import njit, prange
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

STABILITY_OUTFILENAME = "stability_values.normfinder.csv"


############################################################################
# HELPERS
############################################################################


def nanmean(a: np.ndarray, axis: int | None = None) -> np.ndarray:
    """np.nanmean without RuntimeWarnings; all-NaN slices yield NaN."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(a, axis=axis)


############################################################################
# NUMBA-ACCELERATED KERNEL
############################################################################


@njit(parallel=True, cache=True)
def compute_minvars(z: np.ndarray, target_idx: np.ndarray) -> np.ndarray:
    """
    z: (n_genes, n_samples) array, NaN = missing
    target_idx: indices of genes for which to compute the minvar statistic
    returns: 1D array of length len(target_idx)

    For each target gene i, computes min over j != i of the variance of the
    pairwise differences (y_i - y_j) across the samples where BOTH genes are
    valid, divided by 4. Pairs with fewer than 2 shared valid samples are
    skipped; a target gene with no valid sample yields NaN; a target gene with
    no usable pair yields +inf.
    """
    ngenes, nsamples = z.shape

    # should not happen as it is controlled before, but just in case
    if nsamples < 2:
        raise ValueError("Number of samples must be at least 2")

    minvars = np.empty(target_idx.shape[0], dtype=np.float32)
    for k in prange(target_idx.shape[0]):
        i = target_idx[k]
        xi = z[i]

        # checking if counts for this gene are all NaN
        n_valid = 0
        for s in range(nsamples):
            if not np.isnan(xi[s]):
                n_valid += 1
        if n_valid < 1:
            minvars[k] = np.nan
            continue

        minv = np.inf
        for j in range(ngenes):
            if j == i:
                continue
            # single pass over samples: sum and sum of squares of valid diffs
            s1 = 0.0
            s2 = 0.0
            n = 0
            for s in range(nsamples):
                a = xi[s]
                b = z[j, s]
                if not (np.isnan(a) or np.isnan(b)):
                    d = a - b
                    s1 += d
                    s2 += d * d
                    n += 1
            if n < 2:
                continue
            m = s1 / n
            var = (s2 - n * m * m) / (n - 1)
            if var < 0.0:  # guard against catastrophic cancellation
                var = 0.0
            if var < minv:
                minv = var
        minvars[k] = minv / 4.0 if np.isfinite(minv) else np.inf
    return minvars


#####################################################
# NORMFINDER CLASS
#####################################################


@dataclass
class NormFinder:
    count_lf: pl.LazyFrame
    design_df: pl.DataFrame

    genes: list[str] = field(init=False)

    group_to_samples_dict: dict[str, list[str]] = field(init=False)

    n_groups: int = field(init=False)
    n_genes: int = field(init=False)

    # core data (set in __post_init__)
    samples: list[str] = field(init=False)  # column order of X
    matrix: np.ndarray = field(init=False)  # (n_genes, n_samples) float32, NaN = missing

    def __post_init__(self):
        # format_design
        self.design_df = self.design_df.with_columns(
            pl.concat_str([pl.col("batch"), pl.col("condition")], separator="_").alias(
                "group"
            )
        ).select("sample", "group")

        # make dict associating a group to the list of its samples
        group_to_sample_df = self.design_df.group_by("group", maintain_order=True).agg(
            "sample"
        )  # maintain order is better for repeatability and testing

        self.group_to_samples_dict = {
            d["group"]: d["sample"] for d in group_to_sample_df.to_dicts()
        }

        groups = list(self.group_to_samples_dict.keys())
        self.n_groups = len(groups)

        self.genes = (
            self.count_lf.select(config.GENE_ID_COLNAME).collect().to_series().to_list()
        )
        self.n_genes = len(self.genes)

        if self.n_genes <= 2:
            logger.error("Too few genes")
            sys.exit(100)

        # single collect of the count matrix, in design order
        self.samples = self.design_df["sample"].to_list()
        count_df = self.count_lf.select(self.samples).collect()
        self.matrix = (
            count_df
            .cast(pl.Float32)
            .fill_null(np.nan)  # unify nulls and NaNs
            .to_numpy()
        )  # (n_genes, n_samples), row order = self.genes

    def _sample_idx(self, samples: list[str]) -> np.ndarray:
        pos = {s: i for i, s in enumerate(self.samples)}
        return np.array([pos[s] for s in samples], dtype=np.int64)

    def _intragroup_variance(self, Y: np.ndarray) -> np.ndarray:
        """
        Unbiased intragroup variance for one group.

        Y: (n_genes, ng) slice of X for this group's samples.
        Returns (n_genes,) vector.

        r_igj = y_igj - mean(y_ig*) - mean(y_*gj) + mean(y_*g*)
        var_i = sum_j r_igj^2 / ((ng - 1) * (1 - 2 / n_genes))
              - sum_i var_i / (n_genes * (n_genes - 1))
        Genes with a negative estimate are replaced by their minvar / 4.
        """
        ng = Y.shape[1]
        if ng == 1:
            # if only one sample in the group, there's no variance
            return np.zeros(self.n_genes, dtype=np.float64)

        gene_means = nanmean(Y, axis=1)  # mean(y_ig*), (n_genes,)
        sample_means = nanmean(Y, axis=0)  # mean(y_*gj), (ng,)
        overall_mean = nanmean(gene_means)  # mean(y_*g*), scalar

        R = Y - gene_means[:, None] - sample_means[None, :] + overall_mean
        ss = np.nansum(R * R, axis=1)  # sum over samples of r_igj^2
        var = ss / ((ng - 1) * (1.0 - 2.0 / self.n_genes))

        # subtract average over all gene pairs
        total = np.nansum(var)
        var = var - total / (self.n_genes * (self.n_genes - 1))

        # if some values are negative, replace them by minvar / 4
        neg_idx = np.nonzero(np.isfinite(var) & (var < 0.0))[0]
        if neg_idx.size:
            var[neg_idx] = compute_minvars(Y, neg_idx.astype(np.uint32))

        return var

    def _intergroup_variance(
        self, gene_means_in_groups: np.ndarray, dataset_overall_mean: float
    ) -> np.ndarray:
        """
        Unbiased intergroup variance.

        gene_means_in_groups: (n_genes, n_groups) matrix of per-group gene means.
        Returns (n_groups, n_genes) matrix of squared residuals:
        (mean(y_ig*) - mean(y_i**) - mean(y_*g*) + mean(y_***))^2
        """
        mean_over_groups_for_gene = nanmean(gene_means_in_groups, axis=1)  # (n_genes,)
        mean_over_genes_for_group = nanmean(gene_means_in_groups, axis=0)  # (n_groups,)

        R = (
            gene_means_in_groups
            - mean_over_groups_for_gene[:, None]
            - mean_over_genes_for_group[None, :]
            + dataset_overall_mean
        )
        return (R**2).T  # (n_groups, n_genes)

    @staticmethod
    def _gamma_factor(inter: np.ndarray, group_mean_variance: np.ndarray) -> float:
        """Shrinkage factor from the NormFinder model."""
        n_groups, n_genes = inter.shape
        first_term = np.nansum(inter) / ((n_groups - 1) * (n_genes - 1))
        second_term = np.nansum(group_mean_variance) / (n_groups * n_genes)
        return max(first_term - second_term, 0.0)  # set to 0 if negative

    def compute_stability_scoring(self) -> pl.DataFrame:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # UNBIASED INTRAGROUP VARIANCE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logger.info("Computing intragroup variances")
        intra = np.empty((self.n_groups, self.n_genes), dtype=np.float64)
        gene_means_in_groups = np.empty((self.n_genes, self.n_groups), dtype=np.float64)
        group_overall_means = np.empty(self.n_groups, dtype=np.float64)
        n_samples_per_group = np.empty(self.n_groups, dtype=np.float64)

        for g, samples in tqdm(enumerate(self.group_to_samples_dict.values()), total=self.n_groups):
            Y = self.matrix[:, self._sample_idx(samples)]  # (n_genes, ng)
            gene_means = nanmean(Y, axis=1)
            intra[g] = self._intragroup_variance(Y)
            gene_means_in_groups[:, g] = gene_means
            group_overall_means[g] = nanmean(gene_means)
            n_samples_per_group[g] = len(samples)

        logger.info("Adjusting variances by group size")
        group_mean_variance = intra / n_samples_per_group[:, None]  # (n_groups, n_genes)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # INTERGROUP VARIANCE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logger.info("Computing intergroup variances")
        valid_overall_means = group_overall_means[~np.isnan(group_overall_means)]
        dataset_overall_mean = (
            nanmean(valid_overall_means) if valid_overall_means.size else np.nan
        )
        inter = self._intergroup_variance(gene_means_in_groups, dataset_overall_mean)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # STABILITY VALUES
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logger.info("Shrinking intragroup and intergroup variances using gamma factor")
        gamma = self._gamma_factor(inter, group_mean_variance)
        logger.info(f"Gamma factor: {gamma:.6g}")

        difnew = inter * gamma / (gamma + group_mean_variance)
        varnew = group_mean_variance + gamma * group_mean_variance / (
            gamma + group_mean_variance
        )

        logger.info("Computing stability values")
        stability = nanmean(
            np.abs(difnew) + np.sqrt(varnew), axis=0
        )  # mean over groups, (n_genes,)

        return pl.DataFrame(
            {
                config.GENE_ID_COLNAME: self.genes,
                config.NORMFINDER_STABILITY_VALUE_COLNAME: stability,
            }
        )


#####################################################
# FUNCTIONS
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute NormFinder stability values for each gene in the dataset"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--design", type=Path, dest="design_file", required=True, help="Design file"
    )
    return parser.parse_args()


def export_stability(stabilities: pl.DataFrame):
    """Export stability values to CSV file."""
    logger.info(f"Exporting stability values to: {STABILITY_OUTFILENAME}")
    write_csv_with_floats(stabilities, STABILITY_OUTFILENAME, float_precision=5)


def main():
    args = parse_args()

    logger.info(f"Getting counts from {args.count_file}")
    count_lf = pl.scan_parquet(args.count_file)

    logger.info(f"Getting design from {args.design_file}")
    design_df = pl.read_csv(args.design_file)
    # filter design df to keep only samples that are present in the count dataframe
    design_df = design_df.filter(
        pl.col("sample").is_in(count_lf.collect_schema().names())
    )

    nfd = NormFinder(count_lf, design_df)
    stabilities = nfd.compute_stability_scoring()

    logger.info(f"Stability values:\n{stabilities}")
    export_stability(stabilities)


if __name__ == "__main__":
    main()
