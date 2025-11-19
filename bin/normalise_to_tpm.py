#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


TPM_NORM_SUFFIX = ".tpm.csv"


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
        "--annotation",
        type=Path,
        dest="annotation_file",
        required=True,
        help="File containing gene annotations (and gene lengths in particular)",
    )
    return parser.parse_args()


def is_raw_counts(df: pd.DataFrame):
    """Check if the data are raw counts (integers)."""
    return all(df.dtypes.apply(lambda x: pd.api.types.is_integer_dtype(x)))


def is_tpm(df: pd.DataFrame):
    """Check if the data are TPM (sum to 1e6 per sample)."""
    sample_sums = df.sum(axis=0)
    return all((sample_sums - 1e6).abs() < 1e-3)  # Allow for floating-point precision


def is_fpkm_or_rpkm(df: pd.DataFrame):
    """Check if the data are FPKM or RPKM (not raw, not TPM)."""
    return not is_raw_counts(df) and not is_tpm(df)


def process_to_tpm(df: pd.DataFrame, gene_lengths: list):
    """
    Process raw counts, FPKM, or RPKM to TPM.
    - For raw counts: Calculate RPKM, then TPM.
    - For FPKM/RPKM: Convert directly to TPM.
    """
    if is_raw_counts(df):
        # Calculate RPKM
        total_reads = df.sum(axis=0)
        rpkm = df.div(gene_lengths, axis=0) / total_reads * 1e9
        # Convert RPKM to TPM
        tpm = rpkm.div(rpkm.sum(axis=0), axis=1) * 1e6
        return tpm
    elif is_fpkm_or_rpkm(df):
        # Convert FPKM/RPKM to TPM
        tpm = df.div(df.sum(axis=0), axis=1) * 1e6
        return tpm
    elif is_tpm(df):
        print("Data are already TPM. No conversion needed.")
        return df
    else:
        raise ValueError("Could not determine data type.")


def export_count_data(quantile_normalized_counts: pd.DataFrame, count_file: Path):
    """Export gene expression data to CSV files."""
    # replace .csv / .tsv by .tpm.csv
    outfilename = ".".join(count_file.name.split(".")[:-1]) + TPM_NORM_SUFFIX
    logger.info(f"Exporting quantile normalised counts to: {outfilename}")
    quantile_normalized_counts.reset_index().to_parquet(outfilename)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_file = args.count_file

    logger.info(f"Normalising {count_file.name}")
    count_df = pd.read_csv(count_file, index_col=0)
    count_df.index.name = config.ENSEMBL_GENE_ID_COLNAME

    quantile_normalized_counts = quantile_normalise(count_df, args.target_distribution)

    export_count_data(quantile_normalized_counts, count_file)


if __name__ == "__main__":
    main()
