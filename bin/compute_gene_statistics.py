#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from dataclasses import dataclass, field
from pathlib import Path

import config
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
ALL_GENES_RESULT_OUTFILE_SUFFIX = "stats_all_genes.csv"

RCV_MULTIFILER = 1.4826  # see https://pmc.ncbi.nlm.nih.gov/articles/PMC9196089/

# quantile intervals
NB_QUANTILES = 100


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
    count_df: pl.DataFrame
    nb_nulls_per_samples_df: pl.DataFrame
    max_ratio_null_valid_sample: float
    platform: str | None = field(default=None)

    gene_count_per_sample_df: pl.DataFrame = field(init=False)
    stat_df: pl.DataFrame = field(init=False)
    samples: list[str] = field(init=False)
    samples_with_low_gene_count: list[str] = field(init=False)

    def __post_init__(self):
        self.samples = [
            col for col in self.count_df.columns if col != config.GENE_ID_COLNAME
        ]
        self.samples_with_low_gene_count = self.get_samples_with_low_gene_count()

    def get_colname(self, colname: str) -> str:
        return f"{self.platform}_{colname}" if self.platform else colname

    def get_valid_counts(self) -> pl.DataFrame:
        return self.count_df.select(pl.exclude(config.GENE_ID_COLNAME))

    def get_samples_with_low_gene_count(self) -> list[str]:
        return (
            self.nb_nulls_per_samples_df.filter(
                pl.col(config.SAMPLE_COLNAME).is_in(self.samples)
            )
            .filter(pl.col(config.RATIO_COLNAME) > self.max_ratio_null_valid_sample)
            .select(config.SAMPLE_COLNAME)
            .to_series()
            .to_list()
        )

    def get_main_statistics(self) -> pl.DataFrame:
        """
        Compute count descriptive statistics for each gene in the count dataframe.
        """
        logger.info("Getting descriptive statistics")
        # computing main stats
        augmented_count_df = self.count_df.with_columns(
            mean=pl.concat_list(self.samples).row.mean(),
            std=pl.concat_list(self.samples).row.std(),
            median=pl.concat_list(self.samples).row.median(),
            mad=pl.concat_list(self.samples).row.mad(),
        )

        return augmented_count_df.select(
            pl.col(config.GENE_ID_COLNAME),
            pl.col("mean").alias(self.get_colname(config.MEAN_COLNAME)),
            pl.col("std").alias(self.get_colname(config.STANDARD_DEVIATION_COLNAME)),
            pl.col("median").alias(self.get_colname(config.MEDIAN_COLNAME)),
            pl.col("mad").alias(self.get_colname(config.MAD_COLNAME)),
            (pl.col("std") / pl.col("mean")).alias(
                self.get_colname(config.VARIATION_COEFFICIENT_COLNAME)
            ),
            (pl.col("mad") / pl.col("median") * RCV_MULTIFILER).alias(
                self.get_colname(config.ROBUST_COEFFICIENT_OF_VARIATION_MEDIAN_COLNAME)
            ),
        )

    def compute_ratios_null_values(self):
        # the samples showing a low gene count will not be taken into account for the zero count penalty
        valid_samples = [
            sample
            for sample in self.samples
            if sample not in self.samples_with_low_gene_count
        ]

        nb_nulls = self.count_df.select(
            pl.exclude(config.GENE_ID_COLNAME).is_null()
        ).sum_horizontal()

        if valid_samples:
            nb_nulls_valid_samples = self.count_df.select(
                pl.col(valid_samples).is_null()
            ).sum_horizontal()
        else:
            nb_nulls_valid_samples = nb_nulls

        self.stat_df = self.stat_df.with_columns(
            (nb_nulls / len(self.samples)).alias(
                self.get_colname(config.RATIO_NULLS_COLNAME)
            ),
            (nb_nulls_valid_samples / len(valid_samples)).alias(
                self.get_colname(config.RATIO_NULLS_VALID_SAMPLES_COLNAME)
            ),
        )

    def compute_ratio_zeros(self):
        nb_zeros = self.count_df.select(
            pl.exclude(config.GENE_ID_COLNAME) == 0
        ).sum_horizontal()

        self.stat_df = self.stat_df.with_columns(
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
        self.stat_df = self.stat_df.with_columns(
            (pl.col(mean_colname).rank() / pl.col(mean_colname).count() * NB_QUANTILES)
            .floor()
            .cast(pl.Int8)
            # we want the only value = NB_QUANTILES to be NB_QUANTILES - 1
            # because the last quantile interval is [NB_QUANTILES - 1, NB_QUANTILES]
            .replace({NB_QUANTILES: NB_QUANTILES - 1})
            .alias(self.get_colname(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME))
        )

    def compute_statistics(self) -> pl.DataFrame:
        logger.info("Computing statistics and stability score")
        # getting expression statistics
        self.stat_df = self.get_main_statistics()
        # adding column for nb of null values for each gene
        self.compute_ratios_null_values()
        # adding a column for the frequency of zero values
        self.compute_ratio_zeros()
        # getting quantile intervals
        self.get_quantile_intervals()
        return self.stat_df


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
    parser.add_argument(
        "--ratio-nulls-per-sample",
        type=Path,
        dest="ratio_nulls_per_samples",
        required=True,
        help="Ratio of null values per sample",
    )
    parser.add_argument(
        "--max-ratio-null-valid-sample",
        type=float,
        dest="max_ratio_null_valid_sample",
        required=True,
        help="Maximum ratio of null values for a sample to be considered valid",
    )
    parser.add_argument("--platform", type=str, help="Platform name")
    return parser.parse_args()


def get_counts(file: Path) -> pl.DataFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.read_parquet(file).sort(config.GENE_ID_COLNAME, descending=False)


def export_data(stat_df: pl.DataFrame, platform: str | None):
    """Export gene expression data to CSV files."""
    outfile = (
        f"{platform}.{ALL_GENES_RESULT_OUTFILE_SUFFIX}"
        if platform
        else ALL_GENES_RESULT_OUTFILE_SUFFIX
    )
    logger.info(f"Exporting statistics for all genes to: {outfile}")
    stat_df.write_csv(outfile, float_precision=config.CSV_FLOAT_PRECISION)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    # putting all counts into a single dataframe
    logger.info("Loading count data...")
    count_df = get_counts(args.count_file)
    logger.info(
        f"Loaded count data with {count_df.shape[0]} rows and {count_df.shape[1]} columns"
    )

    ratio_nulls_per_samples_df = pl.read_csv(args.ratio_nulls_per_samples)

    # computing statistics (mean, standard deviation, coefficient of variation, quantiles)
    gene_stat = GeneStatistician(
        count_df,
        ratio_nulls_per_samples_df,
        args.max_ratio_null_valid_sample,
        args.platform,
    )
    stat_df = gene_stat.compute_statistics()

    # exporting computed data
    export_data(stat_df, args.platform)


if __name__ == "__main__":
    main()
