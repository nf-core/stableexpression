#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import sys
import polars as pl
from pathlib import Path
import logging
from functools import reduce

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ALL_COUNTS_PARQUET_OUTFILENAME = "all_counts.parquet"

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


def parse_count_file(count_file: Path) -> pl.LazyFrame:
    lf = pl.scan_csv(count_file)
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


def join_dfs(lf1: pl.LazyFrame, lf2: pl.LazyFrame):
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
    return [
        col for col in lf.collect_schema().names() if col != ENSEMBL_GENE_ID_COLNAME
    ]


def get_counts(files: list[Path]) -> pl.LazyFrame:
    """Get all count data from a list of files.

    The files are merged into a single dataframe. The ENSEMBL_GENE_ID_COLNAME column is cast
    to String, and all other columns are cast to Float64.
    """
    # lazy loading
    lfs = get_valid_lazy_dfs(files)
    # joining all count files
    merged_lf = reduce(join_dfs, lfs)

    # checking if filtered count dataframe is empty
    if merged_lf.limit(1).collect().is_empty():
        logger.error("No data found in any of the input count datasets...")
        sys.exit(100)

    # casting count columns to Float64
    # casting gene id column to String
    count_columns = get_count_columns(merged_lf)
    return (
        merged_lf.select(
            [pl.col(ENSEMBL_GENE_ID_COLNAME).cast(pl.String)]
            + [pl.col(column).cast(pl.Float64) for column in count_columns]
        )
        .fill_null(0)
        .fill_nan(0)
    )


def filter_out_genes_not_always_present(count_lf: pl.LazyFrame):
    filtered_count_lf = count_lf.filter(
        pl.concat_list(pl.exclude(ENSEMBL_GENE_ID_COLNAME)).list.min() > 0
    )
    # checking if filtered count dataframe is empty
    if filtered_count_lf.limit(1).collect().is_empty():
        logger.error("No gene left after filtering for expression > 0 in all samples")
        sys.exit(101)

    return filtered_count_lf


def export_count_data(filtered_count_lf: pl.LazyFrame):
    """Export gene expression data to CSV files."""
    logger.info(f"Exporting normalised counts to: {ALL_COUNTS_PARQUET_OUTFILENAME}")
    filtered_count_lf.collect().write_parquet(ALL_COUNTS_PARQUET_OUTFILENAME)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_files = [Path(file) for file in args.count_files.split(" ")]

    # putting all counts into a single dataframe
    count_lf = get_counts(count_files)
    nb_genes = count_lf.select(pl.len()).collect().item()

    filtered_count_lf = filter_out_genes_not_always_present(count_lf)
    nb_genes_left = count_lf.select(pl.len()).collect().item()

    logger.info(
        f"Kept {nb_genes_left} genes ({100 * nb_genes_left / nb_genes:.2f}%) after filtering"
    )

    export_count_data(filtered_count_lf)


if __name__ == "__main__":
    main()
