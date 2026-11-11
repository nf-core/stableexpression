#!/usr/bin/env python3

# NormFinder stability values — scalable NumPy + Numba implementation.
# Corrected to faithfully reproduce the original algorithm:
#   Andersen CL, Jensen JL, Ørntoft TF. Normalization of Real-Time Quantitative
#   Reverse Transcription-PCR Data: A Model-Based Variance Estimation Approach to
#   Identify Genes Suited for Normalization. Cancer Res. 2004;64(15):5245-5250.
#   (Supplementary material, equations B, C and 1.10)
#
# Corrections vs. the previous version:
#   1. The stability value now uses the |d_ig| intergroup terms of equation C
#      (unsquared residuals), not their squares. Gamma is still estimated from
#      squared residuals (as in the model), but the shrunk intergroup component
#      entering the stability value is |d_ig| * gamma/(gamma+sigma^2_ig/n_g).
#      This restores the paper's interpretation: the stability value is an error
#      magnitude on the (log-)expression scale.
#   2. GroupDif (2 * max_g |d_ig|) and GroupSD (pooled intragroup SD) are now
#      computed, matching the official R/Excel outputs.
#   3. Single-group designs are now supported (the paper's "one group only"
#      branch: GroupSD-based stability), instead of dividing by (n_groups - 1).
#   4. The best two-gene combination (equation 1.10 of the supplementary
#      material) is now computed, with the pStabLim candidate pre-filter.
#   5. Overall means are taken over all (gene, group) entries, matching the R
#      reference implementation.
#
# NOTE ON INPUT SCALE: the official NormFinder software log-transforms data by
# default. Feed this script already log-transformed values (e.g. log2 of
# linear-scale quantities), or raw Ct values (which behave like log-scale).
#
# Written originally by Olivier Coen. Corrected and optimised with Mistral.
# Released under the MIT license.

import argparse
import logging
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl
from numba import njit, prange
from tqdm import tqdm

import config
from common import export_parquet

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

STABILITY_OUTFILENAME = "stability_values.normfinder.parquet"
PAIRS_OUTFILENAME = "pair_stability_values.normfinder.parquet"

# Candidate genes with a stability value below this limit are considered for
# the best two-gene combination (pStabLim in the R reference; 0.25 there).
P_STAB_LIM_DEFAULT = 0.25


############################################################################
# HELPERS
############################################################################


def nanmean(values: np.ndarray, axis: int | None = None) -> np.ndarray:
    """np.nanmean without RuntimeWarnings; all-NaN slices yield NaN."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=axis)


def nansum(values: np.ndarray, axis: int | None = None) -> np.ndarray:
    """np.nansum without RuntimeWarnings (all-NaN slices yield 0, as np.nansum)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nansum(values, axis=axis)


def nanmax(values: np.ndarray, axis: int | None = None) -> np.ndarray:
    """np.nanmax without RuntimeWarnings; all-NaN slices yield NaN."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmax(values, axis=axis)


############################################################################
# NUMBA-ACCELERATED KERNEL
############################################################################


@njit(parallel=True, cache=True)
def compute_minvars(
    expression_matrix: np.ndarray,
    target_gene_indices: np.ndarray,
) -> np.ndarray:
    """Compute minvars.

    expression_matrix: (n_genes, n_samples) array, NaN = missing
    target_gene_indices: indices of genes for which to compute the minvar statistic
    returns: 1D array of length len(target_gene_indices)

    For each target gene i, computes min over j != i of the variance of the
    pairwise differences (y_i - y_j) across the samples where BOTH genes are
    valid, divided by 4. Pairs with fewer than 2 shared valid samples are
    skipped; a target gene with no valid sample yields NaN; a target gene with
    no usable pair yields +inf.
    """
    n_genes, n_samples = expression_matrix.shape

    # should not happen as it is controlled before, but just in case
    if n_samples < 2:
        raise ValueError("Number of samples must be at least 2")

    minvars = np.empty(target_gene_indices.shape[0], dtype=np.float32)
    for target_position in prange(target_gene_indices.shape[0]):
        target_gene = target_gene_indices[target_position]
        target_profile = expression_matrix[target_gene]

        # checking if counts for this gene are all NaN
        n_valid = 0
        for sample_index in range(n_samples):
            if not np.isnan(target_profile[sample_index]):
                n_valid += 1
        if n_valid < 1:
            minvars[target_position] = np.nan
            continue

        min_variance = np.inf
        for other_gene in range(n_genes):
            if other_gene == target_gene:
                continue
            # single pass over samples: sum and sum of squares of valid diffs
            diff_sum = 0.0
            squared_diff_sum = 0.0
            n_valid_pairs = 0
            for sample_index in range(n_samples):
                target_value = target_profile[sample_index]
                other_value = expression_matrix[other_gene, sample_index]
                if not (np.isnan(target_value) or np.isnan(other_value)):
                    diff = target_value - other_value
                    diff_sum += diff
                    squared_diff_sum += diff * diff
                    n_valid_pairs += 1
            if n_valid_pairs < 2:
                continue
            diff_mean = diff_sum / n_valid_pairs
            pair_variance = (
                squared_diff_sum - n_valid_pairs * diff_mean * diff_mean
            ) / (n_valid_pairs - 1)
            if pair_variance < 0.0:  # guard against catastrophic cancellation
                pair_variance = 0.0
            min_variance = min(min_variance, pair_variance)
        minvars[target_position] = (
            min_variance / 4.0 if np.isfinite(min_variance) else np.inf
        )
    return minvars


#####################################################
# NORMFINDER CLASS
#####################################################


@dataclass
class NormFinder:
    """Class to compute the NormFinder stability score for each gene in the dataset."""

    count_lf: pl.LazyFrame
    design_df: pl.DataFrame

    genes: list[str] = field(init=False)

    group_to_samples_dict: dict[str, list[str]] = field(init=False)

    n_groups: int = field(init=False)
    n_genes: int = field(init=False)

    # core data (set in __post_init__)
    samples: list[str] = field(init=False)  # column order of expression_matrix
    expression_matrix: np.ndarray = field(init=False)  # (n_genes, n_samples) float32, NaN = missing

    def __post_init__(self):
        # format_design
        self.design_df = self.design_df.with_columns(
            pl.concat_str([pl.col("batch"), pl.col("condition")], separator="_").alias(
                "group",
            ),
        ).select("sample", "group")

        # make dict associating a group to the list of its samples
        group_to_sample_df = self.design_df.group_by("group", maintain_order=True).agg(
            "sample",
        )  # maintain order is better for repeatability and testing

        self.group_to_samples_dict = {
            group_row["group"]: group_row["sample"]
            for group_row in group_to_sample_df.to_dicts()
        }

        groups = list(self.group_to_samples_dict.keys())
        self.n_groups = len(groups)

        self.genes = self.count_lf.select(config.GENE_ID_COLNAME).collect().to_series().to_list()
        self.n_genes = len(self.genes)

        if self.n_genes <= 2:
            logger.error("Too few genes")
            sys.exit(100)

        # single collect of the count matrix, in design order
        self.samples = self.design_df["sample"].to_list()
        count_df = self.count_lf.select(self.samples).collect()
        self.expression_matrix = (
            count_df.cast(pl.Float32)
            .fill_null(np.nan)  # unify nulls and NaNs
            .to_numpy()
        )  # (n_genes, n_samples), row order = self.genes

    def _sample_idx(self, samples: list[str]) -> np.ndarray:
        sample_to_position = {sample: i for i, sample in enumerate(self.samples)}
        return np.array(
            [sample_to_position[sample] for sample in samples], dtype=np.int64,
        )

    def _intragroup_variance(self, group_expression_matrix: np.ndarray) -> np.ndarray:
        """Compute unbiased intragroup variance for one group (equation B).

        group_expression_matrix: (n_genes, n_group_samples) slice of
            expression_matrix for this group's samples.
        Returns (n_genes,) vector.

        r_igj = y_igj - mean(y_ig*) - mean(y_*gj) + mean(y_*g*)
        var_i = sum_j r_igj^2 / ((n_group_samples - 1) * (1 - 2 / n_genes))
              - sum_i var_i / (n_genes * (n_genes - 1))
        Genes with a negative estimate are replaced by their minvar / 4.
        """
        n_group_samples = group_expression_matrix.shape[1]
        if n_group_samples == 1:
            # if only one sample in the group, there's no variance
            return np.zeros(self.n_genes, dtype=np.float64)

        gene_means = nanmean(group_expression_matrix, axis=1)  # mean(y_ig*), (n_genes,)
        sample_means = nanmean(group_expression_matrix, axis=0)  # mean(y_*gj), (n_group_samples,)
        overall_group_mean = nanmean(group_expression_matrix)  # mean(y_*g*), scalar

        centered_expression = (
            group_expression_matrix
            - gene_means[:, None]
            - sample_means[None, :]
            + overall_group_mean
        )
        sum_of_squared_residuals = nansum(
            centered_expression * centered_expression, axis=1,
        )  # sum over samples of r_igj^2
        variance_estimates = sum_of_squared_residuals / (
            (n_group_samples - 1) * (1.0 - 2.0 / self.n_genes)
        )

        # subtract average over all gene pairs
        total_variance = nansum(variance_estimates)
        variance_estimates = variance_estimates - total_variance / (
            self.n_genes * (self.n_genes - 1)
        )

        # if some values are negative, replace them by minvar / 4
        negative_variance_indices = np.nonzero(
            np.isfinite(variance_estimates) & (variance_estimates < 0.0),
        )[0]
        if negative_variance_indices.size:
            variance_estimates[negative_variance_indices] = compute_minvars(
                group_expression_matrix, negative_variance_indices.astype(np.uint32),
            )

        return variance_estimates

    @staticmethod
    def _intergroup_residuals(
        gene_means_in_groups: np.ndarray,
        dataset_overall_mean: float,
    ) -> np.ndarray:
        """Two-way centered intergroup residuals d_ig (equation C, UNSQUARED).

        gene_means_in_groups: (n_genes, n_groups) matrix of per-group gene means.
        Returns (n_genes, n_groups) matrix:
        d_ig = mean(y_ig*) - mean(y_i**) - mean(y_*g*) + mean(y_***)
        """
        mean_over_groups_for_gene = nanmean(gene_means_in_groups, axis=1)  # (n_genes,)
        mean_over_genes_for_group = nanmean(gene_means_in_groups, axis=0)  # (n_groups,)

        return (
            gene_means_in_groups
            - mean_over_groups_for_gene[:, None]
            - mean_over_genes_for_group[None, :]
            + dataset_overall_mean
        )

    @staticmethod
    def _gamma_factor(
        intergroup_residuals: np.ndarray,
        variance_of_group_means: np.ndarray,
    ) -> float:
        """Shrinkage factor gamma from the NormFinder model.

        Estimated from the SQUARED residuals (as in the reference
        implementation), even though the stability value itself uses the
        unsquared |d_ig| terms.
        """
        n_groups, n_genes = intergroup_residuals.T.shape
        if n_groups <= 1 or n_genes <= 1:
            return 0.0
        between_group_term = nansum(intergroup_residuals**2) / (
            (n_groups - 1) * (n_genes - 1)
        )
        within_group_term = nansum(variance_of_group_means) / (n_groups * n_genes)
        return max(between_group_term - within_group_term, 0.0)  # set to 0 if negative

    def _single_group_variance(self, expression_matrix: np.ndarray) -> np.ndarray:
        """Compute  unbiased variance for a single-group design ("one group only").

        expression_matrix: (n_genes, n_samples) matrix, all samples.
        Returns (n_genes,) unbiased variance estimates (sigma^2), with
        negative estimates replaced by minvar / 4.
        """
        n_genes, n_samples = expression_matrix.shape
        if n_samples < 2:
            raise ValueError("Number of samples must be at least 2")

        gene_means = nanmean(expression_matrix, axis=1)  # (n_genes,)
        sample_means = nanmean(expression_matrix, axis=0)  # (n_samples,)
        overall_mean = nanmean(expression_matrix)

        centered_expression = (
            expression_matrix
            - gene_means[:, None]
            - sample_means[None, :]
            + overall_mean
        )
        variance_estimates = nansum(
            centered_expression * centered_expression, axis=1,
        ) / ((n_samples - 1) * (1.0 - 2.0 / n_genes))

        total_variance = nansum(variance_estimates)
        variance_estimates = variance_estimates - total_variance / (
            n_genes * (n_genes - 1)
        )

        # if some values are negative, replace them by minvar / 4
        negative_variance_indices = np.nonzero(
            np.isfinite(variance_estimates) & (variance_estimates < 0.0),
        )[0]
        if negative_variance_indices.size:
            variance_estimates[negative_variance_indices] = compute_minvars(
                expression_matrix, negative_variance_indices.astype(np.uint32),
            )

        return variance_estimates

    def _candidate_genes(self, stability_values: np.ndarray) -> np.ndarray:
        """Compute indices of genes eligible for the two-gene combination analysis.

        Genes with a stability value below p_stab_lim; if fewer than 4
        candidates and there are more than 3 genes, the 4 most stable genes
        are used instead (as in the R reference implementation).
        """
        candidate_indices = np.nonzero(stability_values < P_STAB_LIM_DEFAULT)[0]
        if candidate_indices.size < 4:
            if self.n_genes > 3:
                fourth_smallest_stability = np.sort(stability_values)[3]
                candidate_indices = np.nonzero(
                    stability_values <= fourth_smallest_stability,
                )[0]
            else:
                candidate_indices = np.arange(self.n_genes)
        return candidate_indices

    def compute_stability_scoring(self) -> tuple[pl.DataFrame, list[tuple[str, str, float]]]:
        """Compute stability scores for all genes, and the candidate gene pairs."""
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # UNBIASED INTRAGROUP VARIANCE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logger.info("Computing intragroup variances")
        intragroup_variances = np.empty((self.n_groups, self.n_genes), dtype=np.float64)
        gene_means_in_groups = np.empty((self.n_genes, self.n_groups), dtype=np.float64)
        n_samples_per_group = np.empty(self.n_groups, dtype=np.float64)

        for group_index, group_samples in tqdm(
            enumerate(self.group_to_samples_dict.values()), total=self.n_groups,
        ):
            group_expression_matrix = self.expression_matrix[
                :, self._sample_idx(group_samples),
            ]  # (n_genes, n_group_samples)
            gene_means_in_groups[:, group_index] = nanmean(
                group_expression_matrix, axis=1,
            )
            intragroup_variances[group_index] = self._intragroup_variance(
                group_expression_matrix,
            )
            n_samples_per_group[group_index] = len(group_samples)

        logger.info("Adjusting variances by group size")
        variance_of_group_means = intragroup_variances / n_samples_per_group[
            :, None,
        ]  # (n_groups, n_genes)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # INTERGROUP RESIDUALS (equation C, unsquared)
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logger.info("Computing intergroup residuals")
        # overall mean over all (gene, group) entries, as in the reference
        dataset_overall_mean = float(nanmean(gene_means_in_groups))
        intergroup_residuals = self._intergroup_residuals(
            gene_means_in_groups, dataset_overall_mean,
        )

        if self.n_groups > 1:
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # STABILITY VALUES (multiple groups)
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            logger.info("Shrinking intragroup and intergroup variances using gamma factor")
            gamma = self._gamma_factor(intergroup_residuals, variance_of_group_means)
            logger.info(f"Gamma factor: {gamma:.6g}")

            # shrunk intergroup term: |d_ig| * gamma / (gamma + sigma^2_ig / n_g)
            shrunk_intergroup_terms = (
                np.abs(intergroup_residuals.T)
                * gamma
                / (gamma + variance_of_group_means)
            )
            shrunk_intragroup_variances = variance_of_group_means + (
                gamma * variance_of_group_means
                / (gamma + variance_of_group_means)
            )

            logger.info("Computing stability values")
            stability_values = nanmean(
                shrunk_intergroup_terms + np.sqrt(shrunk_intragroup_variances),
                axis=0,
            )  # mean over groups

            # GroupDif: 2 * max_g |d_ig| (as in the official output)
            group_dif = 2.0 * nanmax(np.abs(intergroup_residuals), axis=1)
            # GroupSD: pooled common intragroup standard deviation
            group_sd = np.sqrt(
                nansum(
                    (n_samples_per_group[:, None] - 1.0) * intragroup_variances,
                    axis=0,
                )
                / (np.sum(n_samples_per_group) - self.n_groups),
            )

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # BEST TWO-GENE COMBINATION (equation 1.10)
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            gene_pairs = self._pair_analysis_multiple_groups(
                intergroup_residuals,
                intragroup_variances,
                n_samples_per_group,
                gamma,
                stability_values,
            )
        else:
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # ONE GROUP ONLY: stability value is the GroupSD
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            logger.info("Single group design: computing GroupSD-based stability values")
            gene_standard_deviations = np.sqrt(
                self._single_group_variance(self.expression_matrix),
            )
            stability_values = gene_standard_deviations
            group_dif = np.zeros(self.n_genes)
            group_sd = gene_standard_deviations

            gene_pairs = self._pair_analysis_single_group(gene_standard_deviations)

        result = pl.DataFrame(
            {
                config.GENE_ID_COLNAME: self.genes,
                config.NORMFINDER_STABILITY_VALUE_COLNAME: stability_values,
                "group_dif": group_dif,
                "group_sd": group_sd,
            },
        ).sort(config.NORMFINDER_STABILITY_VALUE_COLNAME)

        return result, gene_pairs

    def _pair_analysis_multiple_groups(
        self,
        intergroup_residuals: np.ndarray,  # (n_genes, n_groups)
        intragroup_variances: np.ndarray,  # (n_groups, n_genes)
        n_samples_per_group: np.ndarray,  # (n_groups,)
        gamma: float,
        stability_values: np.ndarray,
    ) -> list[tuple[str, str, float]]:
        """Compute combined stability for pairs of candidate genes (eq. 1.10).

        Vectorised port of MakeComb2 from the R reference implementation.
        Requires n_genes > 2 (the combination formula diverges otherwise).
        """
        if self.n_genes <= 2:
            return []

        candidate_indices = self._candidate_genes(stability_values)
        if candidate_indices.size < 2:
            return []

        gene_count_factor = np.sqrt(self.n_genes / (self.n_genes - 2.0))

        gene_pairs = []
        for first_position in range(candidate_indices.size):
            for second_position in range(first_position + 1, candidate_indices.size):
                first_gene = candidate_indices[first_position]
                second_gene = candidate_indices[second_position]
                first_gene_residuals = intergroup_residuals[first_gene]  # (n_groups,)
                second_gene_residuals = intergroup_residuals[second_gene]  # (n_groups,)
                first_gene_variances = intragroup_variances[:, first_gene]  # (n_groups,)
                second_gene_variances = intragroup_variances[:, second_gene]  # (n_groups,)

                # [gamma * d / (gamma + s/n)] terms, averaged over the two genes
                pair_stability_terms = (
                    np.abs(
                        gamma * first_gene_residuals
                        / (gamma + first_gene_variances / n_samples_per_group)
                        + gamma * second_gene_residuals
                        / (gamma + second_gene_variances / n_samples_per_group),
                    )
                    * gene_count_factor
                    / 2.0
                )
                # variance contribution of the averaged pair
                pair_stability_terms = pair_stability_terms + np.sqrt(
                    first_gene_variances / n_samples_per_group
                    + gamma * first_gene_variances
                    / (n_samples_per_group * gamma + first_gene_variances)
                    + second_gene_variances / n_samples_per_group
                    + gamma * second_gene_variances
                    / (n_samples_per_group * gamma + second_gene_variances),
                ) / 2.0
                gene_pairs.append(
                    (
                        self.genes[first_gene],
                        self.genes[second_gene],
                        float(np.mean(pair_stability_terms)),
                    ),
                )

        gene_pairs.sort(key=lambda gene_pair: gene_pair[2])
        return gene_pairs

    def _pair_analysis_single_group(
        self,
        gene_standard_deviations: np.ndarray,
    ) -> list[tuple[str, str, float]]:
        """Compute combined stability for pairs of genes in a single-group design.

        Port of the "one group only" branch of the R reference implementation:
        the two candidate genes are replaced by their per-sample average, and
        the stability of the combined "gene" is recomputed.
        """
        if self.n_genes <= 3:
            return []

        # candidate genes: sd < min(sd) + 0.1 (as in the reference)
        sd_limit = np.nanmin(gene_standard_deviations) + 0.1
        candidate_indices = np.nonzero(gene_standard_deviations < sd_limit)[0]
        if candidate_indices.size < 2:
            return []

        gene_pairs = []
        for first_position in range(candidate_indices.size):
            for second_position in range(first_position + 1, candidate_indices.size):
                first_gene = candidate_indices[first_position]
                second_gene = candidate_indices[second_position]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    combined_profile = np.nanmean(
                        self.expression_matrix[[first_gene, second_gene], :], axis=0,
                    )  # (n_samples,)
                combined_matrix = np.vstack(
                    [
                        np.delete(self.expression_matrix, [first_gene, second_gene], axis=0),
                        combined_profile[None, :],
                    ],
                )  # (n_genes - 1, n_samples)
                combined_gene_sd = np.sqrt(
                    self._single_group_variance(combined_matrix),
                )[-1]  # stability of the combined "gene"
                gene_pairs.append(
                    (
                        self.genes[first_gene],
                        self.genes[second_gene],
                        float(combined_gene_sd),
                    ),
                )

        gene_pairs.sort(key=lambda gene_pair: gene_pair[2])
        return gene_pairs


#####################################################
# FUNCTIONS
#####################################################


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute NormFinder stability values for each gene in the dataset",
    )
    parser.add_argument(
        "--counts",
        type=Path,
        dest="count_file",
        required=True,
        help="Count file (already log-transformed or Ct-scale)",
    )
    parser.add_argument(
        "--design",
        type=Path,
        dest="design_file",
        required=True,
        help="Design file",
    )
    parser.add_argument(
        "--export-pairs",
        action="store_true",
        dest="export_pairs",
        help="Whether to export the two-gene combination stability values",
    )
    return parser.parse_args()


def export_stability(stabilities: pl.DataFrame) -> None:
    """Export stability values to parquet file."""
    logger.info(f"Exporting stability values to: {STABILITY_OUTFILENAME}")
    export_parquet(stabilities, STABILITY_OUTFILENAME)


def export_pairs(gene_pairs: list[tuple[str, str, float]]) -> None:
    """Export the two-gene combination stability values to parquet file."""
    if not gene_pairs:
        logger.info("No eligible two-gene combination to export")
        return
    pairs_df = pl.DataFrame(
        {
            "gene_1": [gene_pair[0] for gene_pair in gene_pairs],
            "gene_2": [gene_pair[1] for gene_pair in gene_pairs],
            "combined_stability": [gene_pair[2] for gene_pair in gene_pairs],
        },
    )
    logger.info(f"Exporting pair stability values to: {PAIRS_OUTFILENAME}")
    export_parquet(pairs_df, PAIRS_OUTFILENAME)


def main():
    args = parse_args()

    logger.info(f"Getting counts from {args.count_file}")
    count_lf = pl.scan_parquet(args.count_file)

    logger.info(f"Getting design from {args.design_file}")
    design_df = pl.read_csv(args.design_file)
    # filter design df to keep only samples that are present in the count dataframe
    design_df = design_df.filter(
        pl.col("sample").is_in(count_lf.collect_schema().names()),
    )

    normfinder = NormFinder(count_lf, design_df)
    stabilities, gene_pairs = normfinder.compute_stability_scoring()

    logger.info(f"Stability values:\n{stabilities}")
    export_stability(stabilities)

    if gene_pairs and args.export_pairs:
        best_pair = gene_pairs[0]
        logger.info(
            f"Best two-gene combination: {best_pair[0]} + {best_pair[1]} "
            f"(combined stability: {best_pair[2]:.4f})",
        )
    export_pairs(gene_pairs)


if __name__ == "__main__":
    main()
