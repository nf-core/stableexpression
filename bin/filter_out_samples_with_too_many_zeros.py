#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import export_parquet, parse_count_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".zeros_filtered.parquet"
RATIO_ZEROS_PER_SAMPLE_OUTFILE = "ratio_zeros_per_sample.csv"
RATIO_ZERO_VALUES_OUTFILE = "ratio_zeros.csv"
NB_REJECTED_SAMPLES_OUTFILE = "nb_rejected_samples.csv"
NB_KEPT_SAMPLES_OUTFILE = "nb_kept_samples.csv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Filter out samples not valid")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--max-zero-ratio",
        type=float,
        dest="max_zero_ratio",
        required=True,
        help="Maximum ratio of zeros allowed",
    )
    return parser.parse_args()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    # putting all counts into a single dataframe
    logger.info("Loading count data...")
    count_df = parse_count_table(args.count_file)
    nb_samples = count_df.shape[1] - 1
    logger.info(
        f"Loaded count data with {len(count_df)} genes and {nb_samples} samples"
    )

    # computing the number of zeros values per sample
    ratio_zeros_df = count_df.select(
        pl.exclude(config.GENE_ID_COLNAME).eq(pl.lit(0)).mean()
    )

    # getting the samples with a zero ratio lower than the max zero ratio
    valid_samples = [
        col
        for col in ratio_zeros_df.columns
        if ratio_zeros_df[col][0] <= args.max_zero_ratio
    ]

    # if at least one valid sample is remaining, making an updated count dataframe
    if valid_samples:
        logger.info(f"Filtered out {count_df.shape[1] - len(valid_samples)} columns")
        valid_count_df = count_df.select(
            pl.col(config.GENE_ID_COLNAME), pl.col(valid_samples)
        )

        outfilename = args.count_file.with_suffix(OUTFILE_SUFFIX).name
        export_parquet(valid_count_df, outfilename)
    else:
        logger.error("No valid columns remaining")

    # collect all ratio values for export
    ratio_values = list(ratio_zeros_df.row(0))
    with open(RATIO_ZERO_VALUES_OUTFILE, "w") as outfile:
        # sorting values in order to having consistent output
        outfile.write(",".join([str(val) for val in sorted(ratio_values)]))

    ratio_zeros_df.write_csv(RATIO_ZEROS_PER_SAMPLE_OUTFILE)

    with open(NB_KEPT_SAMPLES_OUTFILE, "w") as fout:
        fout.write(str(len(valid_samples)))

    with open(NB_REJECTED_SAMPLES_OUTFILE, "w") as fout:
        fout.write(str(nb_samples - len(valid_samples)))

    logger.info("Done")


if __name__ == "__main__":
    main()
