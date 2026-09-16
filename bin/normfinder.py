#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean

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
# POLARS EXTENSIONS
############################################################################


@pl.api.register_expr_namespace("row")
class StatsExtension:
    def __init__(self, expr: pl.Expr):
        self._expr = expr

    def not_null_values(self):
        return self._expr.list.eval(pl.element().drop_nulls().drop_nans()).list

    def mean(self) -> pl.Expr:
        """Mean over non nulls values in row"""
        return self.not_null_values().mean()

    def sum(self) -> pl.Expr:
        """Median over non nulls values in row"""
        return self.not_null_values().sum()

    def min(self) -> pl.Expr:
        """Median over non nulls values in row"""
        return self.not_null_values().min()


############################################################################
# NUMBA-ACCELERATED FUNCTIONS
############################################################################


@njit(parallel=True)
def compute_minvars(z: np.ndarray, target_idx: np.ndarray) -> np.ndarray:
    """
    z: (ngenes, nsamples) array
    target_idx: 1D array of indices (uint32) for which to compute minvar
    returns: 1D array of length len(target_idx)
    """
    ngenes, nsamples = z.shape

    # should not happen as it is controlled before, but just in case
    if nsamples < 2:
        raise ValueError("Number of samples must be at least 2")

    minvars = np.empty(len(target_idx), dtype=np.float32)
    for k in prange(len(target_idx)):
        i = target_idx[k]
        # checking if counts for this gene are all nans
        nb_valid_counts = (~np.isnan(z[i, :])).sum()
        if nb_valid_counts < 1:
            minvars[k] = np.nan
            continue  # skip this gene
        # computing variances of pairwise differences
        minv = 1e18
        for j in prange(ngenes):
            if i == j:
                continue
            diffs = z[i, :] - z[j, :]
            mean = np.sum(diffs) / nsamples  # scalar
            var = np.sum((diffs - mean) ** 2) / (nsamples - 1)  # scalar
            if np.isnan(var):
                continue  # skip
            if var < minv:
                minv = var
        minvars[k] = minv / 4.0 if minv < 1e18 else np.inf
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

    @staticmethod
    def get_overall_mean_for_group(df_with_means_over_samples: pl.DataFrame) -> float:
        return df_with_means_over_samples.mean().item()

    @staticmethod
    def get_means_over_samples(df: pl.DataFrame) -> pl.DataFrame:
        return df.with_columns(
            mean_over_samples_for_gene=pl.concat_list(pl.all()).row.mean()
        ).select("mean_over_samples_for_gene")

    def correct_negative_values(
        self, intra_var_df: pl.DataFrame, group_count_df: pl.DataFrame
    ) -> pl.DataFrame:
        genes_with_negative_values = intra_var_df.select(
            col for col in self.genes if (intra_var_df[col] < 0).all()
        ).columns  # intra_var_df has only one row but it is a dataframe

        # getting indexes of genes for which we must compute minvar
        indexes_of_genes_with_negative_values = np.array(
            [
                i
                for i, gene in enumerate(self.genes)
                if gene in genes_with_negative_values
            ],
            dtype=np.uint32,
        )

        minvars = compute_minvars(
            group_count_df.to_numpy(), indexes_of_genes_with_negative_values
        )

        # associating back minvars to their respective gene
        minvar_dict = {
            gene: minvars[i] for i, gene in enumerate(genes_with_negative_values)
        }
        return intra_var_df.with_columns(
            [pl.lit(val).alias(col) for col, val in minvar_dict.items()]
        )

    def get_unbiased_intragroup_variance_for_group(
        self,
        group_count_df: pl.DataFrame,
        means_over_samples_df: pl.DataFrame,
        group_overall_mean: float,
        samples: list[str],
    ):
        # TODO: see if it's correct
        # if only one sample in the group, there's no variance
        if len(samples) == 1:
            data = {gene: [0] for gene in self.genes}
            return pl.DataFrame(data)

        # lf is a lazyframe with a column being the gene ids (gene_id)
        # and other columns being the samples
        # the current chunk corresponds to only one group
        # means_over_samples_df is a single column dataframe containing the means across each row (ie for each gene across samples)
        ng = len(samples)

        means_over_samples = means_over_samples_df.to_series().rename(
            "mean_over_samples_for_gene"
        )

        mean_over_genes = (
            group_count_df.mean()
            .transpose()
            .to_series()
            .rename("mean_over_genes_for_sample")
        )

        sample_variance_df = (
            group_count_df.hstack(
                [means_over_samples]
            )  # adding column containing means over all samples in this group (for each gene)
            .select(
                [
                    (pl.col(c) - pl.col("mean_over_samples_for_gene")).alias(
                        c
                    )  # y_igj - mean(y_ig*)
                    for c in samples
                ]
            )
            .transpose(
                include_header=True, column_names=self.genes
            )  # columns are now genes
            .hstack(
                [mean_over_genes]
            )  # adding column containing means over all genes (for each sample)
            .select(
                [
                    (
                        (
                            pl.col(c)
                            - pl.col("mean_over_genes_for_sample")
                            + group_overall_mean
                        )
                        ** 2
                    ).alias(
                        c
                    )  # r_igj ^2 = (y_igj - mean(y_ig*) -mean(y_*gj) + mean(y_*g*) ) ^ 2
                    for c in self.genes
                ]
            )
            .transpose(include_header=True, column_names=samples)
            .with_columns(
                sample_variance=pl.concat_list(samples).row.sum()
                / (
                    (ng - 1) * (1 - 2 / self.n_genes)
                )  # sum over j (samples) of r_igj ^2 terms
            )
            .select("sample_variance")
            .transpose()
            .rename({f"column_{i}": gene for i, gene in enumerate(self.genes)})
        )

        # sum of all sample variances for all genes
        sample_variance_sum_over_genes = sample_variance_df.select(
            pl.sum_horizontal(pl.all())
        ).item()  # sum of all s_ij² over all genes

        intra_var_df = sample_variance_df.select(
            [
                (
                    pl.col(c)
                    - sample_variance_sum_over_genes
                    / (self.n_genes * (self.n_genes - 1))
                ).alias(c)
                for c in self.genes
            ]
        )
        # if some values are negative, we need a special process
        corrected_intra_var_df = self.correct_negative_values(
            intra_var_df, group_count_df
        )

        return corrected_intra_var_df

    def get_unbiased_intragroup_variances(self):
        unbiased_intragroup_variance_dfs = []
        means_over_samples_dfs = []
        group_overall_means = []

        for group, samples in tqdm(self.group_to_samples_dict.items()):
            # sub dataframe corresponding to this group
            chunk_df = self.count_lf.select(samples).collect()
            # computing means over samples for each gene
            means_over_samples_df = self.get_means_over_samples(chunk_df)
            # getting overall expression average in the group for all genes
            group_overall_mean = self.get_overall_mean_for_group(means_over_samples_df)

            group_unbiased_intragroup_variance_df = (
                self.get_unbiased_intragroup_variance_for_group(
                    chunk_df, means_over_samples_df, group_overall_mean, samples
                )
            )

            # storing intragroup values for each gene in this group
            unbiased_intragroup_variance_dfs.append(
                group_unbiased_intragroup_variance_df
            )
            # storing means over samples in this group for each gene
            means_over_samples_df = means_over_samples_df.rename(
                {"mean_over_samples_for_gene": group}
            )
            means_over_samples_dfs.append(means_over_samples_df)
            # storing overall mean of expression in this group, for all genes and samples
            group_overall_means.append(group_overall_mean)

        # cast all values to float (to avoid issues when concat)
        unbiased_intragroup_variance_dfs = [
            df.select([pl.col(col).cast(pl.Float32) for col in df.columns])
            for df in unbiased_intragroup_variance_dfs
        ]

        # removing None values in group_overall_means
        # which would originate from group chunk dataframes that are full of null values
        group_overall_means = [mean for mean in group_overall_means if mean is not None]

        # before returning:
        # concatenate together all intragroup variance data to have a single df for all groups
        # stack all means over samples horizontally (becomes a gene * group df )
        # get the mean of group_overall_means to get the overall mean expression value in the count dataframe
        return (
            pl.concat(unbiased_intragroup_variance_dfs),
            pl.concat(means_over_samples_dfs, how="horizontal"),
            mean(group_overall_means),
        )

    def adjust_for_nb_of_samples_in_groups(
        self, unbiased_intragroup_variance_df: pl.DataFrame
    ):
        n_samples_list = [
            len(samples) for samples in self.group_to_samples_dict.values()
        ]
        return unbiased_intragroup_variance_df.with_columns(
            n_samples=pl.Series(n_samples_list)
        ).select([(pl.col(c) / pl.col("n_samples")).cast(pl.Float32).alias(c) for c in self.genes])

    def get_unbiased_intergroup_variance(
        self, gene_means_in_groups_df: pl.DataFrame, dataset_overall_mean: float
    ):
        mean_over_genes = (
            gene_means_in_groups_df.mean()
            .transpose()
            .to_series()
            .rename("mean_over_genes_for_group")
        )

        return (
            gene_means_in_groups_df.with_columns(
                mean_over_groups_for_gene=pl.concat_list(pl.all()).row.mean()
            )
            .select(
                [
                    (pl.col(c) - pl.col("mean_over_groups_for_gene")).alias(c)
                    for c in gene_means_in_groups_df.columns
                ]
            )
            .transpose(column_names=self.genes)
            .hstack([mean_over_genes])
            .select(
                [
                    (
                        pl.col(c)
                        - pl.col("mean_over_genes_for_group")
                        + dataset_overall_mean
                    ).alias(c)
                    for c in self.genes
                ]
            )
            .select(
                [(pl.col(c) ** 2).alias(c) for c in self.genes]
            )  # square to get variance
        )

    def compute_gamma_factor(self, diff_df: pl.DataFrame, vardiff_df: pl.DataFrame):
        logger.info("Computing gamma factor")
        first_term = (
            diff_df.with_columns(
                sum_of_squares=pl.concat_list(pl.all()).row.sum()  # sum over columns
            )
            .select("sum_of_squares")
            .sum()  # sum over rows
            .select(
                (
                    pl.col("sum_of_squares")
                    / ((self.n_groups - 1) * (self.n_genes - 1))
                ).cast(pl.Float32).alias("normalised_sum_of_squares")
            )
            .item()
        )

        second_term = (
            vardiff_df.with_columns(
                sum=pl.concat_list(pl.all()).row.sum()  # sum over columns
            )
            .select("sum")
            .sum()  # sum over rows
            .select(
                (pl.col("sum") / (self.n_groups * self.n_genes)).cast(pl.Float32).alias("normalised_sum")
            )
            .item()
        )

        return max(first_term - second_term, 0)  # set to 0 if negative

    @staticmethod
    def apply_gamma_factor(
        gamma: float, diff_df: pl.DataFrame, vardiff_df: pl.DataFrame
    ):
        difnew = diff_df * gamma / (gamma + vardiff_df)
        varnew = vardiff_df + gamma * vardiff_df / (gamma + vardiff_df)

        return (
            difnew.with_columns(pl.all().cast(pl.Float32)),
            varnew.with_columns(pl.all().cast(pl.Float32))
        )

    def apply_shrinkage(
        self, intergroup_variance_df: pl.DataFrame, group_mean_variance_df: pl.DataFrame
    ):
        gamma = self.compute_gamma_factor(
            intergroup_variance_df, group_mean_variance_df
        )
        return self.apply_gamma_factor(
            gamma, intergroup_variance_df, group_mean_variance_df
        )

    def get_stability_values(
        self, shrunk_intervar_df: pl.DataFrame, shrunk_gr_mean_var_df: pl.DataFrame
    ):
        return (
            (
                shrunk_intervar_df.select([pl.col(c).abs() for c in self.genes])
                + shrunk_gr_mean_var_df.select([pl.col(c).sqrt() for c in self.genes])
            )
            .mean()
            .transpose(
                include_header=True,
                header_name=config.GENE_ID_COLNAME,
                column_names=[config.NORMFINDER_STABILITY_VALUE_COLNAME],
            )
        )

    def compute_stability_scoring(self):
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # UNBIASED INTRAGROUP VARIANCE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logger.info("Computing intragroup variances")
        intragroup_variance_df, gene_means_in_groups_df, dataset_overall_mean = (
            self.get_unbiased_intragroup_variances()
        )

        logger.info("Adjusting variances by group size")
        group_mean_variance_df = self.adjust_for_nb_of_samples_in_groups(
            intragroup_variance_df
        )

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # INTERGROUP VARIANCE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logger.info("Computing intergroup variances")
        intergroup_variance_df = self.get_unbiased_intergroup_variance(
            gene_means_in_groups_df, dataset_overall_mean
        )

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # STABILITY VALUES
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        logger.info("Shrinking intragroup and intergroup variances using gamma factor")
        shrunk_intervar_df, shrunk_gr_mean_var_df = self.apply_shrinkage(
            intergroup_variance_df, group_mean_variance_df
        )

        logger.info("Computing stability values")
        return self.get_stability_values(shrunk_intervar_df, shrunk_gr_mean_var_df)


#####################################################
# FUNCTIONS
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Quantile normalise count data for each sample in the dataset"
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
