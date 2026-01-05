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
ALL_GENE_SUMMARY_OUTFILENAME = "all_genes_summary.csv"
MOST_STABLE_GENE_SUMMARY_OUTFILENAME = "most_stable_genes_summary.csv"
ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME = "all_counts_filtered.parquet"
MOST_STABLE_GENES_COUNTS_OUTFILENAME = (
    "most_stable_genes_transposed_counts_filtered.csv"
)

# nb of top stable genes to select and to display at the end
NB_MOST_STABLE_GENES = 1000
# quantile intervals
NB_QUANTILES = 100
NB_TOP_GENES_TO_SHOW_IN_BOX_PLOTS = 100

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get statistics from count data for each gene"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--stats",
        type=Path,
        dest="stat_file",
        required=True,
        help="File containing statistics for all genes and stability scores by candidate genes",
    )
    parser.add_argument(
        "--platform-stats",
        type=Path,
        dest="platform_stat_files",
        nargs="+",
        help="File containing base statistics for all genes and for all datasets for a specific platform",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        dest="metadata_files",
        help="Metadata file",
    )
    parser.add_argument(
        "--mappings", type=str, dest="mapping_files", help="Mapping file"
    )

    return parser.parse_args()


def parse_stat_file(file: Path) -> pl.DataFrame:
    return pl.read_csv(file).with_columns(
        pl.col(config.GENE_ID_COLNAME).cast(pl.String())
    )


def get_non_empty_dataframes(files: list[Path]) -> list[pl.DataFrame]:
    dfs = [pl.read_csv(file) for file in files]
    return [df for df in dfs if not df.is_empty()]


def cast_cols_to_string(df: pl.DataFrame) -> pl.DataFrame:
    return df.select(
        [pl.col(column).cast(pl.String) for column in df.collect_schema().names()]
    )


def concat_cast_to_string_and_drop_duplicates(files: list[Path]) -> pl.DataFrame:
    """Concatenate DataFrames, cast all columns to String, and drop duplicates.

    The first step is to concatenate the DataFrames. Then, the dataframe is cast
    to String to ensure that all columns have the same data type. Finally, duplicate
    rows are dropped.
    """
    dfs = get_non_empty_dataframes(files)
    dfs = [cast_cols_to_string(df) for df in dfs]
    concat_df = pl.concat(dfs)
    # dropping duplicates
    # casting all columns to String
    return concat_df.unique()


def cast_count_columns_to_float(df: pl.DataFrame) -> pl.DataFrame:
    return df.select(
        pl.col(config.GENE_ID_COLNAME),
        pl.exclude(config.GENE_ID_COLNAME).cast(pl.Float64),
    )


def join_data_on_gene_id(stat_df: pl.DataFrame, *dfs: pl.DataFrame) -> pl.DataFrame:
    """Merge the statistics dataframe with the metadata dataframe and the mapping dataframe."""
    # we need to ensure that the index of stat_df are strings
    for df in dfs:
        stat_df = stat_df.join(df, on=config.GENE_ID_COLNAME, how="left")
    return stat_df


def get_counts(file: Path) -> pl.DataFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.read_parquet(file).sort(config.GENE_ID_COLNAME, descending=False)


def get_metadata(metadata_files: list[Path]) -> pl.DataFrame | None:
    """Retrieve and concatenate metadata from a list of metadata files."""
    if not metadata_files:
        return None
    return concat_cast_to_string_and_drop_duplicates(metadata_files)


def get_mappings(mapping_files: list[Path]) -> pl.DataFrame | None:
    if not mapping_files:
        return None
    concat_df = concat_cast_to_string_and_drop_duplicates(mapping_files)
    # group by new gene IDs and gets the lis
    # convert the list column to a string representation
    # separate the original gene IDs with a semicolon
    return concat_df.group_by(config.GENE_ID_COLNAME).agg(
        pl.col(config.ORIGINAL_GENE_ID_COLNAME)
        .unique()
        .sort()
        .str.join(";")
        .alias(config.ORIGINAL_GENE_IDS_COLNAME)
    )


def get_status(quantile_interval: int) -> str:
    """Return the expression level status of the gene given its quantile interval."""
    if NB_QUANTILES - 5 <= quantile_interval:
        return "Very high expression"
    elif NB_QUANTILES - 10 <= quantile_interval < NB_QUANTILES - 5:
        return "High expression"
    elif 4 < quantile_interval <= 9:
        return "Low expression"
    elif quantile_interval <= 4:
        return "Very low expression"
    else:
        return "Medium range"


def add_expression_level_status(df: pl.DataFrame) -> pl.DataFrame:
    logger.info("Adding expression level status")
    mapping_dict = {
        quantile_interval: get_status(quantile_interval)
        for quantile_interval in range(NB_QUANTILES)
    }
    return df.with_columns(
        pl.col(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
        .replace_strict(mapping_dict)
        .alias(config.EXPRESSION_LEVEL_STATUS_COLNAME)
    )


def get_all_genes_summary(
    stat_summary_df: pl.DataFrame, *dfs: pl.DataFrame
) -> pl.DataFrame:
    """
    Extract the most stable genes from the statistics dataframe.
    """
    # add gene name, description and original gene IDs to statistics summary
    stat_summary_df = join_data_on_gene_id(stat_summary_df, *dfs)
    stat_summary_df = add_expression_level_status(stat_summary_df)
    return stat_summary_df


def get_most_stable_genes_counts(
    log_count_df: pl.DataFrame, stat_summary_df: pl.DataFrame
) -> pl.DataFrame:
    # getting list of top stable genes with their order
    top_genes_with_order = (
        stat_summary_df.head(NB_TOP_GENES_TO_SHOW_IN_BOX_PLOTS)
        .select(config.GENE_ID_COLNAME)
        .with_row_index("sort_order")
    )

    # join to get only existing genes and maintain order
    sorted_transposed_counts_df = log_count_df.join(
        top_genes_with_order, on=config.GENE_ID_COLNAME, how="inner"
    ).sort("sort_order", descending=False)

    # get the actual gene names that were found (in order)
    actual_gene_names = (
        sorted_transposed_counts_df.select(config.GENE_ID_COLNAME).to_series().to_list()
    )
    return sorted_transposed_counts_df.drop(
        ["sort_order", config.GENE_ID_COLNAME]
    ).transpose(column_names=actual_gene_names)


def export_data(
    all_genes_summary_df: pl.DataFrame,
    most_stable_genes_summary_df: pl.DataFrame,
    all_counts_df: pl.DataFrame,
    most_stable_genes_counts_df: pl.DataFrame,
):
    """Export gene expression data to CSV files."""
    logger.info(f"Exporting statistics of all genes to: {ALL_GENE_SUMMARY_OUTFILENAME}")
    all_genes_summary_df.write_csv(
        ALL_GENE_SUMMARY_OUTFILENAME, float_precision=config.CSV_FLOAT_PRECISION
    )

    logger.info(
        f"Exporting statistics of the top stable genes to: {MOST_STABLE_GENE_SUMMARY_OUTFILENAME}"
    )
    most_stable_genes_summary_df.write_csv(
        MOST_STABLE_GENE_SUMMARY_OUTFILENAME, float_precision=config.CSV_FLOAT_PRECISION
    )

    logger.info(f"Exporting all counts to: {ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME}")
    all_counts_df.write_parquet(ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME)

    logger.info(
        f"Exporting counts of the top stable genes to: {MOST_STABLE_GENES_COUNTS_OUTFILENAME}"
    )
    most_stable_genes_counts_df.write_csv(
        MOST_STABLE_GENES_COUNTS_OUTFILENAME, float_precision=config.CSV_FLOAT_PRECISION
    )

    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    metadata_files = (
        [Path(file) for file in args.metadata_files.split(" ")]
        if args.metadata_files is not None
        else []
    )
    mapping_files = (
        [Path(file) for file in args.mapping_files.split(" ")]
        if args.mapping_files is not None
        else []
    )

    count_df = get_counts(args.count_file)

    # getting data, including metadata and mappings
    all_genes_stat_summary_df = parse_stat_file(args.stat_file)

    platform_datasets_stat_dfs = [
        parse_stat_file(file) for file in args.platform_stat_files if file is not None
    ]

    metadata_df = get_metadata(metadata_files)
    mapping_df = get_mappings(mapping_files)
    optional_dfs = [df for df in [metadata_df, mapping_df] if df is not None]

    additional_data_dfs = optional_dfs + platform_datasets_stat_dfs
    all_genes_summary_df = get_all_genes_summary(
        all_genes_stat_summary_df, *additional_data_dfs
    )

    top_stable_stat_summary_df = all_genes_summary_df.head(NB_MOST_STABLE_GENES)

    # reducing dataframe size (it is only used for plotting by MultiQC)
    count_df = cast_count_columns_to_float(count_df)
    most_stable_genes_counts_df = get_most_stable_genes_counts(
        count_df, top_stable_stat_summary_df
    )

    # exporting computed data
    export_data(
        all_genes_summary_df,
        top_stable_stat_summary_df,
        count_df,
        most_stable_genes_counts_df,
    )


if __name__ == "__main__":
    main()
