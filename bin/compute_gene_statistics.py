#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
ALL_GENES_RESULT_OUTFILE_SUFFIX = "stats_all_genes.csv"

RCV_MULTIPLIER = 1.4826  # see https://pmc.ncbi.nlm.nih.gov/articles/PMC9196089/

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
    return parser.parse_args()


def get_counts(file: Path) -> pl.LazyFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.scan_parquet(file).sort(config.GENE_ID_COLNAME, descending=False)


def get_colname(colname: str, platform: str | None) -> str:
    return f"{platform}_{colname}" if platform else colname


def get_samples(lf: pl.LazyFrame) -> list[str]:
    return lf.select(pl.exclude(config.GENE_ID_COLNAME)).collect_schema().names()


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
    lf: pl.LazyFrame, valid_samples: list[str], platform: str | None
) -> pl.LazyFrame:
    
    samples_cols = [col for col in lf.collect_schema().names() if col != config.GENE_ID_COLNAME]
    nb_samples = len(samples_cols) - 1
    found_valid_samples = [sample for sample in valid_samples if sample in samples_cols]
    
    # the samples showing a low gene count will not be taken into account for the zero count penalty
    nb_nulls = (
        lf
        .select(pl.exclude(config.GENE_ID_COLNAME).is_null()) # select all columns except GENE_ID_COLNAME and check if they are null
        .select(pl.sum_horizontal(pl.all()).alias("nb_nulls_all_samples")) # sum the number of null values across all columns
        .collect()
        .to_series()
    )
    
    if found_valid_samples:
        nb_nulls_valid_samples = (
            lf
            .select(pl.col(found_valid_samples).is_null()) # select all columns in valid_samples and check if they are null
            .select(pl.sum_horizontal(pl.all()).alias("nb_nulls_valid_samples")) # sum the number of null values across all columns
            .collect()
            .to_series()
        )
    else:
        nb_nulls_valid_samples = nb_nulls
    
    return lf.select(
        pl.col(config.GENE_ID_COLNAME),
        (nb_nulls / nb_samples).alias(get_colname(config.RATIO_NULLS_COLNAME, platform)),
        (nb_nulls_valid_samples / len(found_valid_samples)).alias(get_colname(config.RATIO_NULLS_VALID_SAMPLES_COLNAME, platform)),
    )


def get_main_statistics(lf: pl.LazyFrame, platform: str | None) -> pl.LazyFrame:
    """
    Compute count descriptive statistics for each gene in the count dataframe.
    """
    logger.info("Getting descriptive statistics")
    samples = get_samples(lf)
    # computing main stats
    augmented_count_lf = lf.with_columns(
        mean=pl.concat_list(samples).row.mean(),
        std=pl.concat_list(samples).row.std(),
        median=pl.concat_list(samples).row.median(),
        mad=pl.concat_list(samples).row.mad(),
    )

    return augmented_count_lf.select(
        pl.col(config.GENE_ID_COLNAME),
        pl.col("mean").alias(get_colname(config.MEAN_COLNAME, platform)),
        pl.col("std").alias(get_colname(config.STANDARD_DEVIATION_COLNAME, platform)),
        pl.col("median").alias(get_colname(config.MEDIAN_COLNAME, platform)),
        pl.col("mad").alias(get_colname(config.MAD_COLNAME, platform)),
        (pl.col("std") / pl.col("mean")).alias(
            get_colname(config.COEFFICIENT_OF_VARIATION_COLNAME, platform)
        ),
        (pl.col("mad") / pl.col("median") * RCV_MULTIPLIER).alias(
            get_colname(config.ROBUST_COEFFICIENT_OF_VARIATION_MEDIAN_COLNAME, platform)
        ),
    )


def compute_ratio_zeros(
    count_lf: pl.LazyFrame, stat_lf: pl.LazyFrame, platform: str
) -> pl.LazyFrame:
    nb_samples = len(get_samples(count_lf))
    nb_zeros_lf = count_lf.select(
        (pl.sum_horizontal(pl.exclude(config.GENE_ID_COLNAME) == 0) / nb_samples).alias(
            get_colname(config.RATIO_ZEROS_COLNAME, platform)
        )
    )
    # return stat_lf
    return pl.concat([stat_lf, nb_zeros_lf], how="horizontal")


def get_quantile_intervals(lf: pl.LazyFrame, platform: str) -> pl.LazyFrame:
    """
    Compute the quantile intervals for the mean expression levels of each gene in the dataframe.

    The function assigns to each gene a quantile interval of its mean cpm compared to all genes.
    """
    logger.info("Getting mean expression quantiles")
    mean_colname = get_colname(config.MEAN_COLNAME, platform)
    return lf.with_columns(
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
        .alias(get_colname(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME, platform))
    )


def export_data(lf: pl.LazyFrame, platform: str | None):
    """Export gene expression data to CSV files."""
    outfile = (
        f"{platform}.{ALL_GENES_RESULT_OUTFILE_SUFFIX}"
        if platform
        else ALL_GENES_RESULT_OUTFILE_SUFFIX
    )
    logger.info(f"Exporting statistics for all genes to: {outfile}")
    lf.sink_csv(outfile, float_precision=config.CSV_FLOAT_PRECISION)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    ratio_nulls_per_samples_df = pl.read_csv(args.ratio_nulls_per_samples)
    valid_samples = get_valid_samples(
        ratio_nulls_per_samples_df, args.max_ratio_null_valid_sample
    )

    logger.info("Loading count data (before missing value imputation)")
    non_imputed_count_lf = get_counts(args.count_file)

    ratio_nulls_lf = compute_ratios_null_values(
        non_imputed_count_lf, valid_samples, args.platform
    )

    # if the user provided an imputed count file, use it; otherwise, use the original count file
    if args.imputed_count_file:
        logger.info("Using imputed count file")
        count_file = args.imputed_count_file
    else:
        logger.info("Using original count file")
        count_file = args.count_file

    logger.info("Loading count data...")
    count_lf = get_counts(count_file)

    logger.info("Computing statistics and stability score")
    # getting expression statistics
    stat_lf = get_main_statistics(count_lf, args.platform)

    # adding column for nb of null values for each gene
    stat_lf = stat_lf.join(
        ratio_nulls_lf, on=config.GENE_ID_COLNAME, how="inner"
    )

    # adding a column for the frequency of zero values
    stat_lf = compute_ratio_zeros(count_lf, stat_lf, args.platform)

    # getting quantile intervals
    stat_lf = get_quantile_intervals(stat_lf, args.platform)

    # exporting computed data
    export_data(stat_lf, args.platform)


if __name__ == "__main__":
    main()
