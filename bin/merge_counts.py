#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from functools import reduce
from pathlib import Path

import config
import polars as pl
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ALL_COUNTS_PARQUET_OUTFILENAME = "all_counts.parquet"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Merge count datasets")
    parser.add_argument(
        "--counts", type=str, dest="count_files", required=True, help="Count files"
    )
    return parser.parse_args()


#####################################################
# COUNTS
#####################################################


def parse_count_file(count_file: Path) -> pl.DataFrame:
    df = pl.read_parquet(count_file)
    # in some cases, the first column may have an empty name or be different than config.GENE_ID_COLNAME
    # in any case, this column must have the config.GENE_ID_COLNAME name
    first_column_name = df.columns[0]
    if first_column_name != config.GENE_ID_COLNAME:
        df = df.rename({first_column_name: config.GENE_ID_COLNAME})
    return df


def is_valid_df(df: pl.DataFrame, file: Path) -> bool:
    """Check if a DataFrame is valid.
    A DataFrame is considered valid if it contains at least one row.
    """
    try:
        return not df.limit(1).is_empty()
    except FileNotFoundError:
        # strangely enough we get this error for some files existing but empty
        logger.error(f"Could not find file {str(file)}")
        return False
    except pl.exceptions.NoDataError as err:
        logger.error(f"File {str(file)} is empty: {err}")
        return False


def get_valid_dfs(files: list[Path]) -> list[pl.DataFrame]:
    """Get a list of valid DataFrames from a list of files.
    A DataFrame is considered valid if it contains at least one row.
    """
    df_dict = {file: parse_count_file(file) for file in tqdm(files)}
    return [df for file, df in df_dict.items()]


def join_count_dfs(df1: pl.DataFrame, df2: pl.DataFrame) -> pl.DataFrame:
    """Join two DataFrames on the config.GENE_ID_COLNAME column.

    The how parameter is set to "full" to include all rows from both dfs.
    The coalesce parameter is set to True to fill NaN values in the
    resulting dataframe with values from the other dataframe.
    """
    return df1.join(df2, on=config.GENE_ID_COLNAME, how="full", coalesce=True)


def get_count_columns(df: pl.DataFrame) -> list[str]:
    """Get all column names except the config.GENE_ID_COLNAME column.

    The config.GENE_ID_COLNAME column contains only gene IDs.
    """
    return df.select(pl.exclude(config.GENE_ID_COLNAME)).columns


def get_counts(files: list[Path]) -> pl.DataFrame:
    """Get all count data from a list of files.

    The files are merged into a single dataframe. The config.GENE_ID_COLNAME column is cast
    to String, and all other columns are cast to Float64.
    """
    logger.info("Parsing counts")
    dfs = get_valid_dfs(files)

    # joining all count files
    logger.info(
        f"Joining count files recursively on the {config.GENE_ID_COLNAME} column"
    )
    merged_df = reduce(join_count_dfs, tqdm(dfs))

    count_columns = get_count_columns(merged_df)
    # casting count columns to Float64
    # casting gene id column to Stringcount_files
    # casting nans to nulls
    logger.info("Cleaning mergeed dataframe")
    return merged_df.select(
        [pl.col(config.GENE_ID_COLNAME).cast(pl.String)]
        + [pl.col(column).cast(pl.Float64) for column in count_columns]
    ).fill_nan(None)


#####################################################
# EXPORT
#####################################################


def export_data(count_df: pl.DataFrame):
    """Export gene expression data."""
    logger.info(f"Exporting normalised counts to: {ALL_COUNTS_PARQUET_OUTFILENAME}")
    count_df.write_parquet(ALL_COUNTS_PARQUET_OUTFILENAME)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_files = [Path(file) for file in args.count_files.split(" ")]
    logger.info(f"Merging {len(count_files)} count files")

    # putting all counts into a single dataframe
    count_df = get_counts(count_files)
    export_data(count_df)


if __name__ == "__main__":
    main()
