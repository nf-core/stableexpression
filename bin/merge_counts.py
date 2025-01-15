#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging
from functools import reduce

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

COUNT_SUMMARY_PARQUET_OUTFILENAME = "all_counts.parquet"

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"


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
    return parser.parse_args()


def is_valid_df(df: pl.LazyFrame, file: Path) -> bool:
    """Check if a LazyFrame is valid.

    A LazyFrame is considered valid if it contains at least one row.
    """
    try:
        return not df.limit(1).collect().is_empty()
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
    df_dict = {file: pl.scan_csv(file) for file in files}
    return [df for file, df in df_dict.items() if is_valid_df(df, file)]


def join_dfs(df1: pl.LazyFrame, df2: pl.LazyFrame):
    """Join two LazyFrames on the ENSEMBL_GENE_ID_COLNAME column.

    The how parameter is set to "full" to include all rows from both dfs.
    The coalesce parameter is set to True to fill NaN values in the
    resulting dataframe with values from the other dataframe.
    """
    return df1.join(df2, on=ENSEMBL_GENE_ID_COLNAME, how="full", coalesce=True)


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


def export_count_data(count_lf: pl.LazyFrame):
    """Export gene expression data to CSV files."""
    logger.info(f"Exporting normalised counts to: {COUNT_SUMMARY_PARQUET_OUTFILENAME}")
    count_lf.collect().write_parquet(COUNT_SUMMARY_PARQUET_OUTFILENAME)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_files = [Path(file) for file in args.count_files.split(" ")]

    # putting all counts into a single dataframe
    count_df = get_counts(count_files)

    export_count_data(count_df)


if __name__ == "__main__":
    main()
