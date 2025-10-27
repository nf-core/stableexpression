#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
from dataclasses import dataclass, field
from typing import ClassVar
import logging

import config

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
ALL_GENES_RESULT_OUTFILE_SUFFIX = "stats_all_genes.csv"


############################################################################
# POLARS EXTENSIONS
############################################################################


@pl.api.register_expr_namespace("row")
class StatsExtension:
    def __init__(self, expr: pl.Expr):
        self._expr = expr

    def not_null_values(self):
        return self._expr.list.drop_nulls().list

    def mean(self) -> pl.Expr:
        """Mean over non nulls values in row"""
        return self.not_null_values().mean()

    def std(self) -> pl.Expr:
        """Std over non nulls values in row"""
        return self.not_null_values().std()

    def median(self) -> pl.Expr:
        """Median over non nulls values in row"""
        return self.not_null_values().median()

    def mad(self) -> pl.Expr:
        """Median Absolute Deviation over non nulls values in row"""
        return (
            self.not_null_values()
            .eval(
                (pl.element() - pl.element().median()).abs().median()
            )  # returns a list with one element
            .list.first()
        )


@dataclass
class GeneStatistician:
    # we want to select samples that show a particularly low nb of genes
    MIN_RATIO_GENE_COUNT_TO_MEAN: ClassVar[float] = 0.75  # experimentally chosen
    # quantile intervals
    NB_QUANTILES: ClassVar[int] = 100

    count_lf: pl.LazyFrame
    platform: str | None = field(default=None)

    gene_count_per_sample_df: pl.DataFrame = field(init=False)
    stat_lf: pl.LazyFrame = field(init=False)
    samples: list[str] = field(init=False)
    samples_with_low_gene_count: list[str] = field(init=False)

    def __post_init__(self):
        self.gene_count_per_sample_df = self.get_gene_counts_per_sample()
        self.samples = (
            self.count_lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME))
            .collect_schema()
            .names()
        )
        self.samples_with_low_gene_count = self.get_samples_with_low_gene_count()

    def get_colname(self, colname: str) -> str:
        return f"{self.platform}_{colname}" if self.platform else colname

    def get_valid_counts(self) -> pl.LazyFrame:
        return self.count_lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME))

    def get_gene_counts_per_sample(self) -> pl.DataFrame:
        """
        Get the number of non-null values per sample.
        :return:
        A polars dataframe containing 2 columns:
            - sample: name of the sample
            - nb_not_nulls: number of non-null values
        """
        return (
            self.count_lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME))
            .count()
            .collect()
            .transpose(
                include_header=True, header_name="sample", column_names=["count"]
            )
        )

    def get_samples_with_low_gene_count(self) -> list[str]:
        mean_gene_count = self.gene_count_per_sample_df[
            config.GENE_COUNT_COLNAME
        ].mean()
        return (
            self.gene_count_per_sample_df.filter(
                (pl.col(config.GENE_COUNT_COLNAME) / mean_gene_count)
                < self.MIN_RATIO_GENE_COUNT_TO_MEAN
            )
            .select(config.SAMPLE_COLNAME)
            .to_series()
            .to_list()
        )

    def get_main_statistics(self) -> pl.LazyFrame:
        """
        Compute count descriptive statistics for each gene in the count dataframe.
        """
        logger.info("Getting descriptive statistics")
        # computing main stats
        augmented_count_lf = self.count_lf.with_columns(
            mean=pl.concat_list(self.samples).row.mean(),
            std=pl.concat_list(self.samples).row.std(),
            median=pl.concat_list(self.samples).row.median(),
            mad=pl.concat_list(self.samples).row.mad(),
        )

        return augmented_count_lf.select(
            pl.col(config.ENSEMBL_GENE_ID_COLNAME),
            pl.col("mean").alias(self.get_colname(config.MEAN_COLNAME)),
            pl.col("std").alias(self.get_colname(config.STANDARD_DEVIATION_COLNAME)),
            pl.col("median").alias(self.get_colname(config.MEDIAN_COLNAME)),
            pl.col("mad").alias(self.get_colname(config.MAD_COLNAME)),
            (pl.col("std") / pl.col("mean")).alias(
                self.get_colname(config.VARIATION_COEFFICIENT_COLNAME)
            ),
        )

    def compute_ratios_null_values(self):
        # the samples showing a low gene count will not be taken into account for the zero count penalty
        valid_samples = [
            sample
            for sample in self.samples
            if sample not in self.samples_with_low_gene_count
        ]

        nb_nulls = (
            self.count_lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME).is_null())
            .collect()
            .sum_horizontal()
        )
        nb_nulls_valid_samples = (
            self.count_lf.select(pl.col(valid_samples).is_null())
            .collect()
            .sum_horizontal()
        )

        self.stat_lf = self.stat_lf.with_columns(
            (nb_nulls / len(self.samples)).alias(
                self.get_colname(config.RATIO_NULLS_COLNAME)
            ),
            (nb_nulls_valid_samples / len(valid_samples)).alias(
                self.get_colname(config.RATIO_NULLS_VALID_SAMPLES_COLNAME)
            ),
        )

    def compute_ratio_zeros(self):
        nb_zeros = (
            self.count_lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME) == 0)
            .collect()
            .sum_horizontal()
        )

        self.stat_lf = self.stat_lf.with_columns(
            (nb_zeros / len(self.samples)).alias(
                self.get_colname(config.RATIO_ZEROS_COLNAME)
            ),
        )

    def get_quantile_intervals(self):
        """
        Compute the quantile intervals for the mean expression levels of each gene in the dataframe.

        The function assigns to each gene a quantile interval of its mean cpm compared to all genes.
        """
        logger.info("Getting cpm quantiles")
        mean_colname = self.get_colname(config.MEAN_COLNAME)
        self.stat_lf = self.stat_lf.with_columns(
            (
                pl.col(mean_colname).rank()
                / pl.col(mean_colname).count()
                * self.NB_QUANTILES
            )
            .floor()
            .cast(pl.Int8)
            # we want the only value = NB_QUANTILES to be NB_QUANTILES - 1
            # because the last quantile interval is [NB_QUANTILES - 1, NB_QUANTILES]
            .replace({self.NB_QUANTILES: self.NB_QUANTILES - 1})
            .alias(self.get_colname(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME))
        )

    def compute_statistics(self) -> pl.LazyFrame:
        logger.info("Computing statistics and stability score")
        # getting expression statistics
        self.stat_lf = self.get_main_statistics()
        # adding column for nb of null values for each gene
        self.compute_ratios_null_values()
        # adding a column for the frequency of zero values
        self.compute_ratio_zeros()
        # getting quantile intervals
        self.get_quantile_intervals()
        return self.stat_lf


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get base statistics from count data for each gene"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument("--platform", type=str, help="Platform name")
    return parser.parse_args()


def get_counts(file: Path) -> pl.LazyFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.scan_parquet(file).sort(config.ENSEMBL_GENE_ID_COLNAME, descending=False)


def export_data(stat_lf: pl.LazyFrame, platform: str):
    """Export gene expression data to CSV files."""
    outfile = f"{platform}.{ALL_GENES_RESULT_OUTFILE_SUFFIX}"
    logger.info(f"Exporting statistics for all genes to: {outfile}")
    stat_lf.collect().write_csv(outfile)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    # putting all counts into a single dataframe
    count_lf = get_counts(args.count_file)

    # computing statistics (mean, standard deviation, coefficient of variation, quantiles)
    gene_stat = GeneStatistician(count_lf, args.platform)
    stat_lf = gene_stat.compute_statistics()

    # exporting computed data
    export_data(stat_lf, args.platform)


if __name__ == "__main__":
    main()
