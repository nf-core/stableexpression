#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
import gzip
from pathlib import Path

import config
import polars as pl
import polars.selectors as cs

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_header(file: Path, sep: str):
    if file.suffix == ".gz":
        fin = gzip.open(file, "rt")
    else:
        fin = open(file, "r")
    header = fin.readline().strip().split(sep)
    first_row = fin.readline().strip().split(sep)

    fin.close()

    if len(header) == len(first_row):
        return header
    elif len(header) == len(first_row) - 1:
        return [config.GENE_ID_COLNAME] + header
    else:
        raise ValueError(
            f"Header has length: {len(header)} while first row has length: {len(first_row)}"
        )


def parse_table(file: Path):
    # parsing header first
    if file.suffix == ".gz":
        ext = file.suffixes[-2]
    else:
        ext = file.suffix
    if ext in [".csv", ".tsv"]:
        # parsing header manually
        sep = "," if ext == ".csv" else "\t"
        header = parse_header(file, sep)
        return pl.read_csv(
            file,
            separator=sep,
            has_header=False,
            skip_rows=1,
            new_columns=header,
            null_values=["NA", "N/A", "na", "n/a"],
        )
    elif ext == ".parquet":
        return pl.read_parquet(file)
    else:
        raise ValueError(f"Unsupported file format: {ext}")


def get_nb_rows(lf: pl.LazyFrame):
    return lf.select(pl.len()).collect().item()


def parse_count_table(file: Path):
    df = parse_table(file)
    first_col = df.columns[0]
    # whatever the name of the first col, rename it to "gene_id"
    return df.rename({first_col: config.GENE_ID_COLNAME}).select(
        pl.col(config.GENE_ID_COLNAME).cast(pl.String()),
        pl.exclude(config.GENE_ID_COLNAME).cast(pl.Float32),
    )


def compute_log2(df: pl.DataFrame) -> pl.DataFrame:
    """
    Compute log2 values.
    """
    return df.select(
        pl.col(config.GENE_ID_COLNAME),
        (pl.exclude(config.GENE_ID_COLNAME) + 1).log(base=2),
    )


def export_parquet(df: pl.DataFrame, outfilename: str):
    logger.info(f"Exporting processed counts to: {outfilename}")
    # round all float columns to avoid inconsistencies during subsequent computations
    # cast float columns to Float32 to fix the
    df.with_columns(cs.float().round(8).cast(pl.Float32)).write_parquet(outfilename)


def write_float_csv(
    df: pl.DataFrame,
    outfilename: str,
    float_precision: int = config.DEFAULT_CSV_FLOAT_PRECISION,
):
    df.write_csv(outfilename, float_precision=float_precision)
