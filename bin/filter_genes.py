#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
"""
Filter out genes with too many zeros across samples
"""

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import export_parquet, parse_count_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".genes_filtered.parquet"
NB_REJECTED_GENES_OUTFILE = "nb_rejected_genes.csv"
NB_KEPT_GENES_OUTFILE = "nb_kept_genes.csv"


#####################################################
# FUNCTIONS
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
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
# MAIN
#####################################################


def main():
    args = parse_args()

    # putting all counts into a single dataframe
    logger.info("Loading count data...")
    count_df = parse_count_table(args.count_file)
    n_genes = count_df.height
    n_samples = count_df.width - 1

    # filtering out genes with too many zeros
    logger.info("Filtering out genes with too many zeros...")
    count_df = count_df.filter(
        pl.sum_horizontal(pl.exclude(config.GENE_ID_COLNAME) == 0) / n_samples <= args.max_zero_ratio
    )

    n_kept_genes = count_df.height
    n_rejected_genes = n_genes - n_kept_genes
    logger.info(f"Kept {n_kept_genes} genes, rejected {n_rejected_genes} genes")

    # exporting
    logger.info("Exporting...")
    outfilename = args.count_file.with_suffix(OUTFILE_SUFFIX).name
    export_parquet(count_df, outfilename)

    with open(NB_KEPT_GENES_OUTFILE, "w") as fout:
        fout.write(str(n_kept_genes))

    with open(NB_REJECTED_GENES_OUTFILE, "w") as fout:
        fout.write(str(n_rejected_genes))

    logger.info("Done")


if __name__ == "__main__":
    main()
