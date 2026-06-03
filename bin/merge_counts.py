#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import hashlib
import json
import logging
from operator import attrgetter
from pathlib import Path

import config
import polars as pl

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


def get_lazyframes(files: list[Path]) -> list[pl.LazyFrame]:
    """Get a list of LazyFrames from a list of files."""
    return [pl.scan_parquet(file, low_memory=True) for file in files]


def get_columns(lf: pl.LazyFrame) -> list[str]:
    return lf.collect_schema().names()


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    return [col for col in get_columns(lf) if col != config.GENE_ID_COLNAME]


def reproducible_hash(lf: pl.LazyFrame) -> str:
    """
    Return a deterministic MD5 hash for a lazyframe.

    Steps:
    1. Convert the tuple (and any nested structures) to a canonical JSON string.
       - `sort_keys=True` guarantees that dictionaries are ordered consistently.
       - `separators=(',', ':')` removes unnecessary whitespace.
    2. Encode the string as UTF‑8 bytes.
    3. Feed the bytes to hashlib.md5 and return the hex digest.

    The result is a 64‑character hexadecimal string that will be identical
    across Python runs, machines, and even different Python versions
    (provided the data types are JSON‑compatible).
    """
    tpl = tuple(get_columns(lf))
    # Canonical JSON representation
    canonical_str = json.dumps(tpl, sort_keys=True, separators=(",", ":"))
    # Encode to bytes
    data_bytes = canonical_str.encode("utf-8")
    # Compute MD5
    hash_obj = hashlib.md5(data_bytes)
    return hash_obj.hexdigest()


def scan_counts(files: list[Path]) -> list[pl.LazyFrame]:
    """
    Get all count data from a list of files.
    """
    logger.info("Parsing counts")
    # sorting them by file name to ensure consistent order between runs
    files.sort(key=attrgetter("name"))

    lfs = get_lazyframes(files)

    # sorting dataframes by a hash on column names
    # this is crucial for consistent output of the script
    # in case multiple files have the same name
    return sorted(lfs, key=lambda lf: reproducible_hash(lf))


def collect_all_gene_ids(lfs: list[pl.LazyFrame]) -> pl.DataFrame:
    """
    Collect all gene IDs from a list of lazyframes.
    """
    logger.info("Getting the full list of gene IDs")
    gene_id_set = set()
    for lf in lfs:
        lf_gene_ids = lf.select(config.GENE_ID_COLNAME).collect().to_series().to_list()
        gene_id_set.update(lf_gene_ids)
    return pl.DataFrame({config.GENE_ID_COLNAME: sorted(list(gene_id_set))})


def make_tmp_sorted_dataframes(
    lfs: list[pl.LazyFrame], gene_id_df: pl.DataFrame
) -> list[Path]:
    """ """
    tmp_files = []
    for i, lf in enumerate(lfs):
        # perform left join from gene ids so that all dataframes can be compared row-wise
        # removing the gene id column for now
        df = gene_id_df.join(
            lf.collect(), on=config.GENE_ID_COLNAME, how="left"
        ).select(pl.exclude(config.GENE_ID_COLNAME))
        outfile = Path(f"tmp.{i}.parquet")
        df.write_parquet(outfile)
        tmp_files.append(outfile)
    return tmp_files


def formating_counts(lf: pl.LazyFrame):
    """
    The config.GENE_ID_COLNAME column is cast
    to String, and all other columns are cast to Float64.
    """
    # casting count columns to Float64
    # casting gene id column to String
    # replacing nans with nulls
    logger.info("Cleaning merged lazyframe")
    return lf.select(
        [pl.col(config.GENE_ID_COLNAME).cast(pl.String)]
        + [pl.col(column).cast(pl.Float64) for column in get_count_columns(lf)]
    ).fill_nan(None)


#####################################################
# EXPORT
#####################################################


def export_data(lf: pl.LazyFrame):
    """Export gene expression data."""
    logger.info(f"Exporting normalised counts to: {ALL_COUNTS_PARQUET_OUTFILENAME}")
    lf.sink_parquet(ALL_COUNTS_PARQUET_OUTFILENAME)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    # parsing count files
    count_files = [Path(file) for file in args.count_files.split(" ")]
    logger.info(f"Merging {len(count_files)} count files")

    lfs = scan_counts(count_files)

    # collecting all gene ids from all lazyframes into a dataframe with one column
    gene_id_df = collect_all_gene_ids(lfs)

    # performing a left join between the sorted list of gene if and each collected lazyframe separately
    # writing this sorted dataframe in a tmp file
    tmp_files = make_tmp_sorted_dataframes(lfs, gene_id_df)

    # scanning the newly created tmp files
    lfs = scan_counts(tmp_files)

    # these files are ready to be merged directly through horizontal concatenation
    # setting strict=True requires all DataFrames to be the same height, raising an error if not.
    merged_lf = pl.concat([gene_id_df.lazy()] + lfs, how="horizontal", strict=True)

    # performing some cleaning / formating operations
    merged_lf = formating_counts(merged_lf)

    # exporting merged data in streaming mode
    export_data(merged_lf)

    # cleaning up tmp files
    for tmp_file in tmp_files:
        tmp_file.unlink()


if __name__ == "__main__":
    main()
