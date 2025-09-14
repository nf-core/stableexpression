#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging
from functools import reduce

import config

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ALL_COUNTS_PARQUET_OUTFILENAME = "all_counts.parquet"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge count datasets"
    )
    parser.add_argument(
        "--counts", type=str, dest="count_files", required=True, help="Count files"
    )
    return parser.parse_args()


#####################################################
# COUNTS
#####################################################


def parse_count_file(count_file: Path) -> pl.LazyFrame:
    lf = pl.scan_parquet(count_file)
    # in some cases, the first column may have an empty name or be different than config.ENSEMBL_GENE_ID_COLNAME
    # in any case, this column must have the config.ENSEMBL_GENE_ID_COLNAME name
    first_column_name = lf.collect_schema().names()[0]
    if first_column_name != config.ENSEMBL_GENE_ID_COLNAME:
        lf = lf.rename({first_column_name: config.ENSEMBL_GENE_ID_COLNAME})
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
    """Join two LazyFrames on the config.ENSEMBL_GENE_ID_COLNAME column.

    The how parameter is set to "full" to include all rows from both dfs.
    The coalesce parameter is set to True to fill NaN values in the
    resulting dataframe with values from the other dataframe.
    """
    return lf1.join(lf2, on=config.ENSEMBL_GENE_ID_COLNAME, how="full", coalesce=True)


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the config.ENSEMBL_GENE_ID_COLNAME column.

    The config.ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME)).collect_schema().names()


def get_counts(files: list[Path]) -> pl.DataFrame:
    """Get all count data from a list of files.

    The files are merged into a single dataframe. The config.ENSEMBL_GENE_ID_COLNAME column is cast
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
            [pl.col(config.ENSEMBL_GENE_ID_COLNAME).cast(pl.String)]
            + [pl.col(column).cast(pl.Float64) for column in count_columns]
        )
        .fill_nan(None)
        .collect()
    )


def get_nb_rows(lf: pl.LazyFrame) -> int:
    return lf.select(pl.len()).collect().item()


#####################################################
# EXPORT
#####################################################


def export_data(count_df: pl.DataFrame ):
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

    # putting all counts into a single dataframe
    count_df = get_counts(count_files)
    export_data(count_df)


if __name__ == "__main__":
    main()
