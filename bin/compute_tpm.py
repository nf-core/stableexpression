#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import polars as pl
from common import compute_log2, parse_count_table, parse_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


OUTFILE_SUFFIX = ".tpm.parquet"

WARNING_REASON_FILE = "warning_reason.txt"
FAILURE_REASON_FILE = "failure_reason.txt"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Normalise data to TPM")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--gene-lengths",
        type=Path,
        dest="gene_lengths_file",
        required=True,
        help="Gene lengths file (CSV format)",
    )
    return parser.parse_args()


def try_cast_to_int(df: pl.DataFrame) -> pl.DataFrame:
    """Try casting columns to integers."""
    count_columns = df.select(pl.exclude(config.GENE_ID_COLNAME)).columns
    # try casting to handle integer values that are float-formated like 1.0
    for col in count_columns:
        is_all_integers = df.select(pl.col(col).round().eq(pl.col(col)).all()).item()
        if is_all_integers:
            df = df.with_columns(pl.col(col).cast(pl.Int64()))
    return df


def is_raw_counts(df: pl.DataFrame) -> bool:
    """Check if the data are raw counts (integers)."""
    count_columns = df.select(pl.exclude(config.GENE_ID_COLNAME)).columns
    return all(
        dtype
        in (
            pl.Int8(),
            pl.Int16(),
            pl.Int32(),
            pl.Int64(),
            pl.UInt8(),
            pl.UInt16(),
            pl.UInt32(),
            pl.UInt64(),
        )
        for dtype in df.select(count_columns).schema.values()
    )


def is_tpm(df: pl.DataFrame) -> bool:
    """Check if the data are TPM (sum to 1e6 per sample)."""
    sample_sums_df = df.select(pl.exclude(config.GENE_ID_COLNAME).sum())
    # a small error is possible, and we assume that if the sum is close to 1e6, it is TPM
    # setting the tolerance to 100
    is_tpm_col_df = sample_sums_df.select((pl.all() - 1e6).abs() < 1e2)
    return is_tpm_col_df.select(
        pl.any_horizontal(pl.all())
    ).item()  # Allow for floating-point precision


def compute_rpkm(df: pl.DataFrame, cdna_length_df: pl.DataFrame) -> pl.DataFrame:
    """
    Process raw counts to RPKM.
    """
    logger.info("Computing RPKM.")
    df = df.join(cdna_length_df, on=config.GENE_ID_COLNAME)
    return df.select(
        pl.col(config.GENE_ID_COLNAME),
        pl.exclude([config.GENE_ID_COLNAME, config.CDNA_LENGTH_COLNAME]).truediv(
            pl.col(config.CDNA_LENGTH_COLNAME)
        ),
    )


def compute_tpm_from_rpkm(rpkm_df: pl.DataFrame) -> pl.DataFrame:
    """
    Process RPKM to TPM.
    """
    logger.info("Computing TPM from RPKM.")
    sums = rpkm_df.select(pl.exclude(config.GENE_ID_COLNAME).sum())
    # Divide each column by its sum and multiply by 1e6
    count_columns = rpkm_df.select(pl.exclude(config.GENE_ID_COLNAME)).columns
    return rpkm_df.select(
        [pl.col(config.GENE_ID_COLNAME)]
        + [(pl.col(col) / sums[col][0] * 1e6).alias(col) for col in count_columns],
    )


def compute_tpm(df: pl.DataFrame, cdna_length_df: pl.DataFrame) -> pl.DataFrame:
    """
    Process raw counts, FPKM, or RPKM to TPM.
    """
    if is_raw_counts(df):
        logger.info("Raw counts detected → computing TPM directly.")
        rpkm_df = compute_rpkm(df, cdna_length_df)
        return compute_tpm_from_rpkm(rpkm_df)
    elif is_tpm(df):
        logger.info("Data are already TPM. No conversion needed.")
        return df
    else:
        # Convert FPKM/RPKM to TPM
        logger.info("Assuming FPKM/RPKM normalisation.")
        return compute_tpm_from_rpkm(df)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    try:
        logger.info("Parsing data")
        count_df = parse_count_table(args.count_file)
        cdna_length_df = parse_table(args.gene_lengths_file)

        logger.info("Converting data types")
        count_df = try_cast_to_int(count_df)

        logger.info(f"Normalising {args.count_file.name}")
        count_df = compute_tpm(count_df, cdna_length_df)

        logger.info("Computing log2 values")
        count_df = compute_log2(count_df)

        outfilename = args.count_file.with_suffix(OUTFILE_SUFFIX).name
        logger.info(f"Exporting TPM normalised counts to: {outfilename}")
        count_df.write_parquet(outfilename)

    except Exception as e:
        logger.error(f"Error occurred while normalising data: {e}")
        msg = "UNEXPECTED ERROR"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)


if __name__ == "__main__":
    main()
