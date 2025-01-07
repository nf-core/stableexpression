#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging
from functools import reduce

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ORIGINAL_GENE_ID_COLNAME = "original_gene_id"
ORIGINAL_GENE_IDS_COLNAME = "original_gene_ids"
ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
MOST_STABLE_GENES_RESULT_OUTFILENAME = "stats_most_stable_genes.csv"
ALL_GENES_RESULT_OUTFILENAME = "stats_all_genes.csv"
COUNT_SUMMARY_OUTFILENAME = "count_summary.csv"
GENE_NAME_COLNAME = "name"
GENE_DESCRIPTION_COLNAME = "description"

VARIATION_COEFFICIENT_COLNAME = "variation_coefficient"
STANDARD_DEVIATION_COLNAME = "standard_deviation"
MEAN_COLNAME = "mean"
EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME = "expression_level_quantile_interval"
QUANTILE_INTERVAL_STATUS_COLNAME = "quantile_interval_status"

NB_QUANTILES = 20
MAX_SELECTED_STABLE_GENES_PER_QUANTILE_INTERVAL = 10

MOST_STABLE_GENES_RESULT_COLS = [
    ENSEMBL_GENE_ID_COLNAME,
    STANDARD_DEVIATION_COLNAME,
    VARIATION_COEFFICIENT_COLNAME,
    MEAN_COLNAME,
    EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME,
    QUANTILE_INTERVAL_STATUS_COLNAME,
    GENE_NAME_COLNAME,
    GENE_DESCRIPTION_COLNAME,
    ORIGINAL_GENE_IDS_COLNAME,
]

ALL_GENES_STATS_COLS = [
    ENSEMBL_GENE_ID_COLNAME,
    MEAN_COLNAME,
    STANDARD_DEVIATION_COLNAME,
    VARIATION_COEFFICIENT_COLNAME,
]

# Get variation coefficient from count data for each gene
# The variation coefficient is the ratio of the standard deviation to the mean
# We want genes that are neither expressed too much nor too little
# Metadata (name and description) are used to annotate the results
# Likewise, mappings (original gene ids) are used to better associate gene ids with their original gene ids
# Usage:
# python get_variation_coefficients.py --counts <count_file> --metadata <metadata_file> --mapping <mapping_file>


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
        "--counts", type=str, dest="count_files", required=True, help="Count file"
    )
    parser.add_argument(
        "--metadata",
        type=str,
        dest="metadata_files",
        required=True,
        help="Metadata file",
    )
    parser.add_argument(
        "--mappings", type=str, dest="mapping_files", required=True, help="Mapping file"
    )
    return parser.parse_args()


def is_valid_df(df: pl.LazyFrame, file: Path) -> bool:
    """Check if a LazyFrame is valid.

    A LazyFrame is considered valid if it contains at least one row.
    """
    try:
        return not df.limit(1).collect().is_empty()
    except (
        FileNotFoundError
    ):  # strangely enough we get this error for files existing but empty
        logger.error(f"Could not find file {str(file)}")
        return False
    except pl.exceptions.NoDataError as err:
        logger.error(f"File {str(file)} is empty: {err}")
        return False


def get_valid_lazy_dfs(files: list[Path]) -> list[pl.LazyFrame]:
    """Get a list of valid LazyFrames from a list of files.

    A LazyFrame is considered valid if it contains at least one row.
    """
    df_dict = {file: pl.scan_csv(file) for file in files}
    return [df for file, df in df_dict.items() if is_valid_df(df, file)]


def join_dfs(df1: pl.LazyFrame, df2: pl.LazyFrame):
    """Join two LazyFrames on the ENSEMBL_GENE_ID_COLNAME column.

    The how parameter is set to "full" to include all rows from both dfs.
    The coalesce parameter is set to True to fill NaN values in the
    resulting dataframe with values from the other dataframe.
    """
    return df1.join(df2, on=ENSEMBL_GENE_ID_COLNAME, how="full", coalesce=True)


def concat_cast_to_string_and_drop_duplicates(files: list[Path]) -> pl.LazyFrame:
    """Concatenate LazyFrames, cast all columns to String, and drop duplicates.

    The first step is to concatenate the LazyFrames. Then, the dataframe is cast
    to String to ensure that all columns have the same data type. Finally, duplicate
    rows are dropped.
    """
    dfs = get_valid_lazy_dfs(files)
    concat_df = pl.concat(dfs)
    # dropping duplicates
    # casting all columns to String
    return concat_df.unique().with_columns(
        [
            pl.col(column).cast(pl.String)
            for column in concat_df.collect_schema().names()
        ]
    )


def get_count_columns(df: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col for col in df.collect_schema().names() if col != ENSEMBL_GENE_ID_COLNAME
    ]


def get_counts(files: list[Path]) -> pl.LazyFrame:
    """Get all count data from a list of files.

    The files are merged into a single dataframe. The ENSEMBL_GENE_ID_COLNAME column is cast
    to String, and all other columns are cast to Float64.
    """
    # lazy loading
    dfs = get_valid_lazy_dfs(files)
    # joining all count files
    merged_df = reduce(join_dfs, dfs)
    # casting count columns to Float64
    # casting gene id column to String
    count_columns = get_count_columns(merged_df)
    return merged_df.with_columns(
        [pl.col(column).cast(pl.Float64) for column in count_columns]
    ).with_columns([pl.col(ENSEMBL_GENE_ID_COLNAME).cast(pl.String)])


def get_metadata(metadata_files: list[Path]) -> pl.LazyFrame:
    """Retrieve and concatenate metadata from a list of metadata files."""
    return concat_cast_to_string_and_drop_duplicates(metadata_files)


def get_mappings(mapping_files: list[Path]) -> pl.LazyFrame:
    concat_df = concat_cast_to_string_and_drop_duplicates(mapping_files)
    # group by new gene IDs and gets the lis
    """Group by new gene IDs, get the list of distinct original gene IDs and convert to a string representation."""
    # t of distinct original gene IDs for each group
    # convert the list column to a string representation
    # separate the original gene IDs with a semicolon
    return concat_df.group_by(ENSEMBL_GENE_ID_COLNAME).agg(
        pl.col(ORIGINAL_GENE_ID_COLNAME)
        .unique()
        .sort()
        .str.join(";")
        .alias(ORIGINAL_GENE_IDS_COLNAME)
    )


def preprocess_count(count_df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Preprocess the count data by replacing NaN values with 0,
    then computing the log2 of the counts.
    """
    # handling NA values (genes that are not found in all datasets)
    # replace NaN values with 0
    count_df = count_df.fill_null(0)

    # getting log2 counts
    count_df = transform_counts_to_log_counts(count_df)

    # sorting dataframe (necessary to get consistent output)
    return count_df.sort(ENSEMBL_GENE_ID_COLNAME, descending=False)


def transform_counts_to_log_counts(count_df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Transform count data to log2 scale.

    The function converts the count data to log2 scale by adding 1 to each count
    (to avoid log(0) and stabilize variance), and then taking the log base 2 of the result.
    The ENSEMBL_GENE_ID_COLNAME column is not transformed, and the transformed columns
    retain their original names.
    """
    logger.info("Getting average log2 counts")
    count_columns = get_count_columns(count_df)
    transformed_cols = [pl.col(ENSEMBL_GENE_ID_COLNAME)] + [
        (pl.col(column) + 1).log(2).alias(column) for column in count_columns
    ]
    return count_df.select(transformed_cols)


def get_main_statistics(count_df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute log2 count descriptive statistics for each gene in the count dataframe.
    """
    logger.info("Getting descriptive statistics")
    count_columns = get_count_columns(count_df)
    return count_df.with_columns(
        mean=pl.concat_list(count_columns).list.mean(),
        std=pl.concat_list(count_columns).list.std(),
    ).select(
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        pl.col("mean").alias(MEAN_COLNAME),
        pl.col("std").alias(STANDARD_DEVIATION_COLNAME),
        (pl.col("std") / pl.col("mean")).alias(VARIATION_COEFFICIENT_COLNAME),
    )


def get_quantile_intervals(df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute the quantile intervals for the mean expression levels of each gene in the dataframe.

    The function assigns to each gene a quantile interval of its mean cpm compared to all genes.
    """
    logger.info("Getting average log2 cpm quantiles")
    return df.with_columns(
        (pl.col(MEAN_COLNAME).rank() / pl.col(MEAN_COLNAME).count() * NB_QUANTILES)
        .floor()
        .cast(pl.Int8)
        # we want the only value = NB_QUANTILES to be NB_QUANTILES - 1
        # because the last quantile interval is [NB_QUANTILES - 1, NB_QUANTILES]
        .replace({NB_QUANTILES: NB_QUANTILES - 1})
        .alias(EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
    )


def compute_statistics(count_df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute descriptive statistics and quantile intervals for gene expression data.

    This function calculates the main statistics (mean, standard deviation, and
    variation coefficient) for each gene. It then computes the quantile intervals based on mean expression levels,
    assigning each gene a quantile interval of its mean counts per million (cpm) compared to all genes.
    """
    # getting expression statistics
    stat_df = get_main_statistics(count_df)

    # getting quantile intervals
    stat_df = get_quantile_intervals(stat_df)

    return stat_df


def merge_data(
    stat_df: pl.LazyFrame, metadata_df: pl.LazyFrame, mapping_df: pl.LazyFrame
) -> pl.LazyFrame:
    """Merge the statistics dataframe with the metadata dataframe and the mapping dataframe."""
    # we need to ensure that the index of stat_df are strings
    return (
        stat_df.join(metadata_df, on=ENSEMBL_GENE_ID_COLNAME, how="left")
        .join(mapping_df, on=ENSEMBL_GENE_ID_COLNAME, how="left")
        .unique()  # just in case
        .sort(STANDARD_DEVIATION_COLNAME, descending=False)  # VERY IMPORTANT
    )


def get_status(quantile_interval: int) -> str:
    """Return the expression level status of the gene given its quantile interval."""
    if quantile_interval == NB_QUANTILES - 1:
        return "very_high_expression"
    elif quantile_interval == NB_QUANTILES - 2:
        return "high_expression"
    elif quantile_interval == 0:
        return "very_low_expression"
    elif quantile_interval == 1:
        return "low_expression"
    else:
        return "ok"


def get_most_stable_genes(stat_df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Extract the most stable genes from the statistics dataframe.

    This function groups the genes by their expression level quantile intervals
    and selects the top most stable genes from each quantile interval.
    The number of genes selected per interval is determined by the constant
    `MAX_SELECTED_STABLE_GENES_PER_QUANTILE_INTERVAL`.
    """
    logger.info("Getting most stable genes per quantile interval")
    return (
        stat_df.group_by(  # the df should be already sorted
            EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME
        )
        .head(MAX_SELECTED_STABLE_GENES_PER_QUANTILE_INTERVAL)
        .with_columns(
            pl.col(EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
            .map_elements(get_status)
            .alias(QUANTILE_INTERVAL_STATUS_COLNAME)
        )
        .sort(
            [
                EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME,
                STANDARD_DEVIATION_COLNAME,
                ENSEMBL_GENE_ID_COLNAME,  # if genes have same quantile and standard deviation, we sort them anyway by gene ID
            ],
            descending=[True, False, False],
        )
        .select([column for column in MOST_STABLE_GENES_RESULT_COLS])
    )


def format_all_genes_dataframe(stat_df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Format the dataframe containing statistics for all genes by selecting the right columns
    and sorting the dataframe by gene ID.
    """
    return stat_df.select([column for column in ALL_GENES_STATS_COLS]).sort(
        ENSEMBL_GENE_ID_COLNAME, descending=False
    )


def export_data(
    most_stable_genes_df: pl.LazyFrame,
    all_genes_stat_df: pl.LazyFrame,
    count_df: pl.LazyFrame,
):
    """Export gene expression data to CSV files."""
    logger.info(
        f"Exporting statistics for the most stable genes to: {MOST_STABLE_GENES_RESULT_OUTFILENAME}"
    )
    most_stable_genes_df.collect().write_csv(MOST_STABLE_GENES_RESULT_OUTFILENAME)

    logger.info(
        f"Exporting statistics for all genes to: {ALL_GENES_RESULT_OUTFILENAME}"
    )
    all_genes_stat_df.collect().write_csv(ALL_GENES_RESULT_OUTFILENAME)

    logger.info(f"Exporting normalised counts to: {COUNT_SUMMARY_OUTFILENAME}")
    count_df.collect().write_csv(COUNT_SUMMARY_OUTFILENAME)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_files = [Path(file) for file in args.count_files.split(" ")]
    metadata_files = [Path(file) for file in args.metadata_files.split(" ")]
    mapping_files = [Path(file) for file in args.mapping_files.split(" ")]

    # putting all counts into a single dataframe
    count_df = get_counts(count_files)
    # preprocessing counts (removing 0 counts / log transformation)
    count_df = preprocess_count(count_df)

    # getting metadata and mappings
    metadata_df = get_metadata(metadata_files)
    mapping_df = get_mappings(mapping_files)

    # computing statistics (mean, standard deviation, coefficient of variation, quantiles)
    stat_df = compute_statistics(count_df)

    stat_df = merge_data(stat_df, metadata_df, mapping_df)

    most_stable_genes_df = get_most_stable_genes(stat_df)

    all_genes_stat_df = format_all_genes_dataframe(stat_df)
    export_data(most_stable_genes_df, all_genes_stat_df, count_df)


if __name__ == "__main__":
    main()
