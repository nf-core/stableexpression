#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging
from functools import reduce

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ALL_COUNTS_PARQUET_OUTFILENAME = "all_counts.parquet"
GENE_COUNT_STATS_OUTFILENAME = "gene_count_statistics.csv"
SKEWNESS_STATS_OUTFILENAME = "skewness_statistics.csv"
KS_TEST_STATS_OUTFILENAME = "ks_test_statistics.csv"
CANDIDATE_GENE_COUNTS_PARQUET_OUTFILENAME = "candidate_gene_counts.parquet"
DISTRIBUTION_CORRELATIONS_OUTFILENAME = "distribution_correlations.csv"

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
STATISTIC_TYPE_COLNAME = "stat_type"
GENE_COUNT_COLNAME = "count"
SKEWNESS_COLNAME = "skewness"
KS_TEST_COLNAME = "kolmogorov_smirnov_pvalue"
SAMPLE_COLNAME = "sample"

STAT_COLNAME_TO_PARAMS = {
    GENE_COUNT_COLNAME: {
        "outfilename": GENE_COUNT_STATS_OUTFILENAME,
        "descending": False,
    },
    SKEWNESS_COLNAME: {"outfilename": SKEWNESS_STATS_OUTFILENAME, "descending": False},
    KS_TEST_COLNAME: {"outfilename": KS_TEST_STATS_OUTFILENAME, "descending": True},
}


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get variation from count data for each gene"
    )
    parser.add_argument(
        "--counts", type=str, dest="count_files", required=True, help="Count files"
    )
    parser.add_argument(
        "--stats",
        type=str,
        dest="dataset_stat_files",
        required=True,
        help="Dataset stats files",
    )
    parser.add_argument(
        "--nb-candidate-genes",
        type=int,
        dest="nb_candidate_genes",
        required=True,
        help="Number of candidate genes to keep",
    )
    return parser.parse_args()


#####################################################
# COUNTS
#####################################################


def parse_count_file(count_file: Path) -> pl.LazyFrame:
    lf = pl.scan_parquet(count_file)
    # in some cases, the first column may have an empty name or be different than ENSEMBL_GENE_ID_COLNAME
    # in any case, this column must have the ENSEMBL_GENE_ID_COLNAME name
    first_column_name = lf.collect_schema().names()[0]
    if first_column_name != ENSEMBL_GENE_ID_COLNAME:
        lf = lf.rename({first_column_name: ENSEMBL_GENE_ID_COLNAME})
    return lf


def is_valid_df(lf: pl.LazyFrame, file: Path) -> bool:
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


def get_valid_lazy_dfs(files: list[Path]) -> list[pl.LazyFrame]:
    """Get a list of valid LazyFrames from a list of files.

    A LazyFrame is considered valid if it contains at least one row.
    """
    lf_dict = {file: parse_count_file(file) for file in files}
    return [lf for file, lf in lf_dict.items() if is_valid_df(lf, file)]


def join_count_dfs(lf1: pl.LazyFrame, lf2: pl.LazyFrame) -> pl.LazyFrame:
    """Join two LazyFrames on the ENSEMBL_GENE_ID_COLNAME column.

    The how parameter is set to "full" to include all rows from both dfs.
    The coalesce parameter is set to True to fill NaN values in the
    resulting dataframe with values from the other dataframe.
    """
    return lf1.join(lf2, on=ENSEMBL_GENE_ID_COLNAME, how="full", coalesce=True)


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME)).collect_schema().names()


def get_counts(files: list[Path]) -> pl.DataFrame:
    """Get all count data from a list of files.

    The files are merged into a single dataframe. The ENSEMBL_GENE_ID_COLNAME column is cast
    to String, and all other columns are cast to Float64.
    """
    # lazy loading
    lfs = get_valid_lazy_dfs(files)
    # joining all count files
    merged_lf = reduce(join_count_dfs, lfs)

    count_columns = get_count_columns(merged_lf)
    # casting count columns to Float64
    # casting gene id column to String
    # casting nans to nulls
    return (
        merged_lf.select(
            [pl.col(ENSEMBL_GENE_ID_COLNAME).cast(pl.String)]
            + [pl.col(column).cast(pl.Float64) for column in count_columns]
        )
        .fill_nan(None)
        .collect()
    )


def get_nb_rows(lf: pl.LazyFrame) -> int:
    return lf.select(pl.len()).collect().item()


#####################################################
# STATISTICS
#####################################################


def parse_stat_file(stat_file: Path) -> pl.DataFrame:
    return pl.read_csv(stat_file, has_header=True)


def merge_stats(stat_files: list[Path]) -> pl.DataFrame:
    stat_dfs = [parse_stat_file(stat_file) for stat_file in stat_files]
    return pl.concat(stat_dfs, how="vertical")


def compute_distances_to_mean(count_df: pl.DataFrame) -> pl.DataFrame:
    corr_dict = {"sample": [], "correlation": []}

    count_df = count_df.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME))
    mean_series = count_df.mean_horizontal()

    for sample in count_df.columns:
        correlation = count_df.select(pl.corr(count_df[sample], mean_series))
        corr_dict["sample"].append(sample)
        corr_dict["correlation"].append(correlation.item())

    return (
        pl.DataFrame(corr_dict)
        .fill_nan(None)
        .sort(by="correlation", descending=True, nulls_last=True)
    )


#####################################################
# CANDIDATE GENES
#####################################################


def get_candidate_gene_counts(
    count_df: pl.DataFrame, nb_candidate_genes: int
) -> pl.DataFrame:
    candidate_gene_lf = (
        count_df.with_columns(
            std=pl.concat_list(pl.exclude(ENSEMBL_GENE_ID_COLNAME))
            .list.drop_nulls()
            .list.std()
        )
        .sort("std", descending=False)
        .head(nb_candidate_genes)
    )
    candidate_gene_ids = (
        candidate_gene_lf.select(ENSEMBL_GENE_ID_COLNAME).to_series().to_list()
    )
    return count_df.filter(pl.col(ENSEMBL_GENE_ID_COLNAME).is_in(candidate_gene_ids))


#####################################################
# EXPORT
#####################################################


def export_data(
    count_df: pl.DataFrame,
    candidate_gene_counts_df: pl.DataFrame,
    corr_df: pl.DataFrame,
):
    """Export gene expression data."""
    logger.info(f"Exporting normalised counts to: {ALL_COUNTS_PARQUET_OUTFILENAME}")
    count_df.write_parquet(ALL_COUNTS_PARQUET_OUTFILENAME)

    logger.info(
        f"Exporting candidate gene counts to: {CANDIDATE_GENE_COUNTS_PARQUET_OUTFILENAME}"
    )
    candidate_gene_counts_df.write_parquet(CANDIDATE_GENE_COUNTS_PARQUET_OUTFILENAME)

    logger.info(
        f"Exporting distribution correlations to: {DISTRIBUTION_CORRELATIONS_OUTFILENAME}"
    )
    corr_df.write_csv(DISTRIBUTION_CORRELATIONS_OUTFILENAME, include_header=False)


def export_individual_statistics(dataset_stats_df: pl.DataFrame):
    for data_col, params in STAT_COLNAME_TO_PARAMS.items():
        outfilename = params["outfilename"]
        logger.info(f"Exporting {data_col} statistics to: {outfilename}")
        sorted_data = dataset_stats_df[[SAMPLE_COLNAME, data_col]].sort(
            data_col, descending=params["descending"]
        )
        sorted_data.write_csv(outfilename, include_header=False)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_files = [Path(file) for file in args.count_files.split(" ")]
    dataset_stat_files = [Path(file) for file in args.dataset_stat_files.split(" ")]

    # putting all counts into a single dataframe
    count_df = get_counts(count_files)
    # putting all stats data into a single dataframe
    dataset_stats_df = merge_stats(dataset_stat_files)

    candidate_gene_counts_df = get_candidate_gene_counts(
        count_df, args.nb_candidate_genes
    )

    # adding stat about divergence to mean distribution
    corr_df = compute_distances_to_mean(count_df)

    export_data(count_df, candidate_gene_counts_df, corr_df)
    export_individual_statistics(dataset_stats_df)


if __name__ == "__main__":
    main()
