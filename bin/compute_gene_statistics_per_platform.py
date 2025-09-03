#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import sys
import polars as pl
from pathlib import Path
from dataclasses import dataclass, field
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# nb of top stable genes to select and to display at the end
DEFAULT_NB_TOP_STABLE_GENES = 1000
# we want to select samples that show a particularly low nb of genes
MIN_RATIO_GENE_COUNT_TO_MEAN = 0.75  # experimentally chosen
WEIGHT_RATIO_NB_NULLS = 1


# outfile names
ALL_GENES_RESULT_OUTFILENAME = "stats_all_genes.csv"

# column names
ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
VARIATION_COEFFICIENT_COLNAME = "variation_coefficient"
STANDARD_DEVIATION_COLNAME = "standard_deviation"
MEAN_COLNAME = "mean"
EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME = "expression_level_quantile_interval"
EXPRESSION_LEVEL_STATUS_COLNAME = "expression_level_status"
GENE_COUNT_COLNAME = "count"
SAMPLE_COLNAME = "sample"
NB_NULLS_COLNAME = "total_nb_nulls"
NB_NULLS_VALID_SAMPLES_COLNAME = "nb_nulls_valid_samples"
NB_ZEROS_COLNAME = "nb_zeros"
STABILITY_SCORE_COLNAME = "stability_score"


# quantile intervals
NB_QUANTILES = 100

NB_TOP_GENES_TO_SHOW_IN_LOG_COUNTS = 100


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get base statistics from count data for each gene. Excludes aberrant datasets."
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--platform", type=str, required=True, help="Platform name"
    )
    return parser.parse_args()


def is_valid_lf(lf: pl.LazyFrame, file: Path) -> bool:
    """Check if a LazyFrame is valid.

    A LazyFrame is considered valid if it contains at least one row.
    """
    try:
        return not lf.limit(1).collect().is_empty()
    except FileNotFoundError:
        # strangely enough we get this error for some files existing but empty
        logger.error(f"Could not find file {str(file)}")
        return False
    except pl.exceptions.NoDataError as err:
        logger.error(f"File {str(file)} is empty: {err}")
        return False


def get_valid_lazy_lfs(files: list[Path]) -> list[pl.LazyFrame]:
    """Get a list of valid LazyFrames from a list of files.

    A LazyFrame is considered valid if it contains at least one row.
    """
    lf_dict = {file: pl.scan_csv(file) for file in files}
    return [lf for file, lf in lf_dict.items() if is_valid_lf(lf, file)]


def cast_cols_to_string(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        [pl.col(column).cast(pl.String) for column in lf.collect_schema().names()]
    )


def concat_cast_to_string_and_drop_duplicates(files: list[Path]) -> pl.LazyFrame:
    """Concatenate LazyFrames, cast all columns to String, and drop duplicates.

    The first step is to concatenate the LazyFrames. Then, the dataframe is cast
    to String to ensure that all columns have the same data type. Finally, duplicate
    rows are dropped.
    """
    lfs = get_valid_lazy_lfs(files)
    lfs = [cast_cols_to_string(lf) for lf in lfs]
    concat_lf = pl.concat(lfs)
    # dropping duplicates
    # casting all columns to String
    return concat_lf.unique()


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME)).collect_schema().names()


def cast_count_columns_to_float32(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        [pl.col(ENSEMBL_GENE_ID_COLNAME)]
        + [pl.col(column).cast(pl.Float32) for column in get_count_columns(lf)]
    )


def get_counts(file: Path) -> pl.LazyFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.scan_parquet(file).sort(ENSEMBL_GENE_ID_COLNAME, descending=False)


def merge_data(
    stat_lf: pl.LazyFrame, metadata_lf: pl.LazyFrame, mapping_lf: pl.LazyFrame
) -> pl.LazyFrame:
    """Merge the statistics dataframe with the metadata dataframe and the mapping dataframe."""
    # we need to ensure that the index of stat_lf are strings
    return stat_lf.join(metadata_lf, on=ENSEMBL_GENE_ID_COLNAME, how="left").join(
        mapping_lf, on=ENSEMBL_GENE_ID_COLNAME, how="left"
    )



def export_data(
    stat_lf: pl.LazyFrame, platform: str
):
    """Export gene expression data to CSV files."""
    outfilename = f"{platform}_{ALL_GENES_RESULT_OUTFILENAME}"
    logger.info(
        f"Exporting statistics for all genes to: {outfilename}"
    )
    stat_lf.collect().write_csv(outfilename)
    logger.info("Done")


#####################################################
#####################################################
# CLASSES
#####################################################
#####################################################


@dataclass
class StabilityScorer:

    count_lf: pl.LazyFrame
    platform: str

    mean_colname: str = field(init=False)
    std_colname: str = field(init=False)
    var_coeff_colname: str = field(init=False)
    nb_nulls_colname: str = field(init=False)
    nb_nulls_valid_samples_colname: str = field(init=False)

    gene_count_per_sample_df: pl.DataFrame = field(init=False)
    stat_lf: pl.LazyFrame = field(init=False)
    count_columns: list[str] = field(init=False)
    samples_with_low_gene_count: list[str] = field(init=False)
    exp_level_quantile_interval_colname: str = field(init=False)
    stability_score_colname: str = field(init=False)

    def __post_init__(self):
        self.count_columns = get_count_columns(self.count_lf)
        self.gene_count_per_sample_df = self.get_gene_counts_per_sample()
        self.samples_with_low_gene_count = self.get_samples_with_low_gene_count()

        self.mean_colname = f"{self.platform}_{MEAN_COLNAME}"
        self.std_colname = f"{self.platform}_{STANDARD_DEVIATION_COLNAME}"
        self.var_coeff_colname = f"{self.platform}_{VARIATION_COEFFICIENT_COLNAME}"
        self.nb_nulls_colname = f"{self.platform}_{NB_NULLS_COLNAME}"
        self.nb_nulls_valid_samples_colname = f"{self.platform}_{NB_NULLS_VALID_SAMPLES_COLNAME}"
        self.exp_level_quantile_interval_colname = f"{self.platform}_{EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME}"
        self.stability_score_colname = f"{self.platform}_{STABILITY_SCORE_COLNAME}"

    def get_valid_counts(self) -> pl.LazyFrame:
        return self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME))

    def get_gene_counts_per_sample(self) -> pl.DataFrame:
        """
        Get the number of non-null values per sample.
        :return:
        A polars dataframe containing 2 columns:
            - sample: name of the sample
            - nb_not_nulls: number of non-null values
        """
        return (
            self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME))
            .count()
            .collect()
            .transpose(
                include_header=True, header_name="sample", column_names=["count"]
            )
        )

    def get_samples_with_low_gene_count(self) -> list[str]:
        mean_gene_count = self.gene_count_per_sample_df[GENE_COUNT_COLNAME].mean()
        return (
            self.gene_count_per_sample_df.filter(
                (pl.col(GENE_COUNT_COLNAME) / mean_gene_count)
                < MIN_RATIO_GENE_COUNT_TO_MEAN
            )
            .select(SAMPLE_COLNAME)
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
            mean=pl.concat_list(self.count_columns).list.drop_nulls().list.mean(),
            std=pl.concat_list(self.count_columns).list.drop_nulls().list.std(),
        )
        return augmented_count_lf.select(
            pl.col(ENSEMBL_GENE_ID_COLNAME),
            pl.col("mean").alias(self.mean_colname),
            pl.col("std").alias(self.std_colname),
            (pl.col("std") / pl.col("mean")).alias(self.var_coeff_colname),
        )

    def compute_nb_null_values(self):
        # the samples showing a low gene count will not be taken into account for the zero count penalty
        cols_to_exclude = [ENSEMBL_GENE_ID_COLNAME] + self.samples_with_low_gene_count
        total_nb_nulls = (
            self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME).is_null())
            .collect()
            .sum_horizontal()
        )
        nb_nulls_valid_samples = (
            self.count_lf.select(pl.exclude(cols_to_exclude).is_null())
            .collect()
            .sum_horizontal()
        )
        self.stat_lf = self.stat_lf.with_columns(
            total_nb_nulls.alias(self.nb_nulls_colname),
            nb_nulls_valid_samples.alias(self.nb_nulls_valid_samples_colname),
        )

    def get_quantile_intervals(self):
        """
        Compute the quantile intervals for the mean expression levels of each gene in the dataframe.

        The function assigns to each gene a quantile interval of its mean cpm compared to all genes.
        """
        logger.info("Getting cpm quantiles")
        self.stat_lf = self.stat_lf.with_columns(
            (pl.col(self.mean_colname).rank() / pl.col(self.mean_colname).count() * NB_QUANTILES)
            .floor()
            .cast(pl.Int8)
            # we want the only value = NB_QUANTILES to be NB_QUANTILES - 1
            # because the last quantile interval is [NB_QUANTILES - 1, NB_QUANTILES]
            .replace({NB_QUANTILES: NB_QUANTILES - 1})
            .alias(self.exp_level_quantile_interval_colname)
        )

    def compute_stability_score(self):
        logger.info("Computing stability score")
        nb_valid_samples = self.gene_count_per_sample_df.select(pl.len()).item() - len(
            self.samples_with_low_gene_count
        )
        ratio_nb_nulls = (
            self.stat_lf.select(
                pl.col(self.nb_nulls_valid_samples_colname) / nb_valid_samples
            )
            .collect()
            .to_series()
        )
        expr = (
            pl.col(self.std_colname) + ratio_nb_nulls * WEIGHT_RATIO_NB_NULLS
        )
        self.stat_lf = self.stat_lf.with_columns(expr.alias(self.stability_score_colname))

    def compute_statistics_and_score(self) -> pl.LazyFrame:
        logger.info("Computing statistics and stability score")
        # getting expression statistics
        self.stat_lf = self.get_main_statistics()
        # adding column for nb of null values for each gene
        self.compute_nb_null_values()
        # computing stability score
        self.compute_stability_score()
        # getting quantile intervals
        self.get_quantile_intervals()

        return self.stat_lf


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
    stability_scorer = StabilityScorer(count_lf, args.platform)
    stat_lf = stability_scorer.compute_statistics_and_score()

    # exporting computed data
    export_data( stat_lf, args.platform )


if __name__ == "__main__":
    main()
