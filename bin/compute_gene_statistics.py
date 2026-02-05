#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from dataclasses import dataclass, field
from pathlib import Path

import config
import polars as pl
from common import write_float_csv
from resource_management import set_max_memory

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
        "--imputed-counts",
        type=Path,
        dest="imputed_count_file",
        help="Count file with imputed missing values",
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
    parser.add_argument(
        "--memory", type=str, dest="memory", required=True, help="Memory in GB"
    )
    return parser.parse_args()


def get_counts(file: Path) -> pl.DataFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.read_parquet(file).sort(config.GENE_ID_COLNAME, descending=False)


def get_colname(colname: str, platform: str | None) -> str:
    return f"{platform}_{colname}" if platform else colname


def get_valid_samples(
    ratio_nulls_per_samples_df: pl.DataFrame, max_ratio_null_valid_sample: float
) -> list[str]:
    """
    Get samples whose ratio of null values is below the maximum ratio.
    """
    return (
        ratio_nulls_per_samples_df.filter(
            pl.col(config.RATIO_COLNAME) <= max_ratio_null_valid_sample
        )
        .select(config.SAMPLE_COLNAME)
        .to_series()
        .to_list()
    )


def compute_ratios_null_values(
    df: pl.DataFrame, valid_samples: list[str], platform: str | None
):
    # the samples showing a low gene count will not be taken into account for the zero count penalty
    nb_nulls = df.select(pl.exclude(config.GENE_ID_COLNAME).is_null()).sum_horizontal()

    if valid_samples:
        nb_nulls_valid_samples = df.select(
            pl.col(valid_samples).is_null()
        ).sum_horizontal()
    else:
        nb_nulls_valid_samples = nb_nulls

    nb_samples = len(df.columns) - 1
    return df.select(
        pl.col(config.GENE_ID_COLNAME),
        (nb_nulls / nb_samples).alias(
            get_colname(config.RATIO_NULLS_COLNAME, platform)
        ),
        (nb_nulls_valid_samples / len(valid_samples)).alias(
            get_colname(config.RATIO_NULLS_VALID_SAMPLES_COLNAME, platform)
        ),
    )


def export_data(stat_df: pl.DataFrame, platform: str | None):
    """Export gene expression data to CSV files."""
    outfile = (
        f"{platform}.{ALL_GENES_RESULT_OUTFILE_SUFFIX}"
        if platform
        else ALL_GENES_RESULT_OUTFILE_SUFFIX
    )
    logger.info(f"Exporting statistics for all genes to: {outfile}")
    write_float_csv(stat_df, outfile)
    logger.info("Done")


#####################################################
#####################################################
# GeneStatistician CLASS
#####################################################
#####################################################


@dataclass
class GeneStatistician:
    count_df: pl.DataFrame
    ratio_nulls_df: pl.DataFrame
    platform: str | None = field(default=None)

    stat_df: pl.DataFrame = field(init=False)
    samples: list[str] = field(init=False)

    def __post_init__(self):
        self.samples = self.count_df.select(pl.exclude(config.GENE_ID_COLNAME)).columns

    def get_colname(self, colname: str) -> str:
        return get_colname(colname, self.platform)

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
                self.get_colname(config.COEFFICIENT_OF_VARIATION_COLNAME)
            ),
            (pl.col("mad") / pl.col("median") * RCV_MULTIFILER).alias(
                self.get_colname(config.ROBUST_COEFFICIENT_OF_VARIATION_MEDIAN_COLNAME)
            ),
        )

    def add_ratio_null_values(self):
        self.stat_df = self.stat_df.join(
            self.ratio_nulls_df, on=config.GENE_ID_COLNAME, how="inner"
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
        logger.info("Getting mean expression quantiles")
        mean_colname = self.get_colname(config.MEAN_COLNAME)
        self.stat_df = self.stat_df.with_columns(
            (
                pl.col(mean_colname).rank(method="ordinal")
                / pl.col(mean_colname).count()
                * NB_QUANTILES
            )
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
        self.add_ratio_null_values()
        # adding a column for the frequency of zero values
        self.compute_ratio_zeros()
        # getting quantile intervals
        self.get_quantile_intervals()
        return self.stat_df


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    set_max_memory(args.memory)

    ratio_nulls_per_samples_df = pl.read_csv(args.ratio_nulls_per_samples)
    valid_samples = get_valid_samples(
        ratio_nulls_per_samples_df, args.max_ratio_null_valid_sample
    )

    logger.info("Loading count data (before missing value imputation)")
    non_imputed_count_df = get_counts(args.count_file)

    ratio_nulls_df = compute_ratios_null_values(
        non_imputed_count_df, valid_samples, args.platform
    )

    # deleting non_imputed_count_df in order to free unused memory
    del non_imputed_count_df

    # if the user provided an imputed count file, use it; otherwise, use the original count file
    if args.imputed_count_file:
        logger.info("Using imputed count file")
        imputed_count_file = args.imputed_count_file
    else:
        logger.info("Using original count file")
        imputed_count_file = args.count_file

    logger.info("Loading count data...")
    count_df = get_counts(imputed_count_file)
    logger.info(
        f"Loaded count data with {count_df.shape[0]} rows and {count_df.shape[1]} columns"
    )

    # computing statistics (mean, standard deviation, coefficient of variation, quantiles)
    gene_stat = GeneStatistician(
        count_df,
        ratio_nulls_df,
        args.platform,
    )
    stat_df = gene_stat.compute_statistics()

    # exporting computed data
    export_data(stat_df, args.platform)


if __name__ == "__main__":
    main()
