#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


TPM_NORM_SUFFIX = ".tpm.parquet"


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
        "--lengths",
        type=Path,
        dest="lengths_file",
        required=True,
        help="File containing gene lengths",
    )
    return parser.parse_args()


def is_raw_counts(df: pd.DataFrame):
    """Check if the data are raw counts (integers)."""
    return all(df.dtypes.apply(lambda x: pd.api.types.is_integer_dtype(x)))


def is_tpm(df: pd.DataFrame):
    """Check if the data are TPM (sum to 1e6 per sample)."""
    sample_sums = df.sum(axis=0)
    return all((sample_sums - 1e6).abs() < 1e-2)  # Allow for floating-point precision


def is_fpkm_or_rpkm(df: pd.DataFrame):
    """Check if the data are FPKM or RPKM (not raw, not TPM)."""
    return not is_raw_counts(df) and not is_tpm(df)


def compute_tpm(df: pd.DataFrame, cdna_length_series: pd.Series):
    """
    Process raw counts, FPKM, or RPKM to TPM.
    """
    if is_raw_counts(df):
        logger.info("Raw counts detected → computing TPM directly.")
        rpk = df.div(cdna_length_series, axis=0)  # read per kilobase
        tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
        return tpm
    elif is_fpkm_or_rpkm(df):
        # Convert FPKM/RPKM to TPM
        logger.info("FPKM/RPKM detected → computing TPM.")
        tpm = df.div(df.sum(axis=0), axis=1) * 1e6
        return tpm
    elif is_tpm(df):
        logger.info("Data are already TPM. No conversion needed.")
        return df
    else:
        raise ValueError("Could not determine data type.")


def parse_data(
    count_file: Path, gene_length_file: Path
) -> tuple[pd.DataFrame, pd.Series]:
    count_df = pd.read_csv(count_file, index_col=0)
    count_df.index.name = config.GENE_ID_COLNAME

    cdna_length_df = pd.read_csv(
        gene_length_file,
        names=[config.GENE_ID_COLNAME, config.CDNA_LENGTH_COLNAME],
    )
    # merge with gene length and extracts it afterwards
    # so that the genes are in the same order in df and in cdna_length_series
    count_df = pd.merge(
        count_df, cdna_length_df, how="left", left_index=True, right_on="gene_id"
    )
    cdna_length_series = count_df[config.CDNA_LENGTH_COLNAME].astype(float)
    count_df = count_df.drop(
        columns=[config.CDNA_LENGTH_COLNAME]
    )  # keep only expression values
    return count_df, cdna_length_series


def export_normalised_data(count_df: pd.DataFrame, count_file: Path):
    """Export gene expression data to Parquet."""
    # replace .csv / .tsv by .tpm.csv
    outfilename = ".".join(count_file.name.split(".")[:-1]) + TPM_NORM_SUFFIX
    logger.info(f"Exporting TPM normalised counts to: {outfilename}")
    count_df.reset_index().to_parquet(outfilename)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing data")
    count_df, cdna_length_series = parse_data(args.count_file, args.gene_length_file)

    logger.info(f"Normalising {args.count_file.name}")

    count_df = compute_tpm(count_df, cdna_length_series)

    export_normalised_data(count_df, args.count_file)


if __name__ == "__main__":
    main()
