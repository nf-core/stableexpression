#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import polars as pl
from common import parse_count_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

CLEANED_COUNTS_SUFFIX = ".cleaned.parquet"
CLEANED_GENE_IDS_SUFFIX = ".cleaned_gene_ids.txt"

FAILURE_REASON_FILE = "failure_reason.txt"

##################################################################
# FUNCTIONS
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Rename gene IDs using mapped IDs")
    parser.add_argument(
        "--count-file", type=Path, required=True, help="Input file containing counts"
    )
    return parser.parse_args()


def clean_ensembl_gene_id_versioning(df: pl.DataFrame):
    """
    Clean Ensembl gene IDs by removing version numbers.
    Remove the dot and the numbers after it in IDs like ENSG00000000003.17
    """
    return df.with_columns(
        pl.when(pl.col(config.GENE_ID_COLNAME).str.starts_with("ENSG"))
        .then(pl.col(config.GENE_ID_COLNAME).str.extract(r"^(ENSG[a-zA-Z0-9]+)", 1))
        .otherwise(pl.col(config.GENE_ID_COLNAME))
        .alias(config.GENE_ID_COLNAME)
    )


def clean_mirna_ids(df: pl.DataFrame):
    """
    Clean miRNA IDs by removing the 5p / 3p identifier.
    """
    return df.with_columns(
        pl.when(pl.col(config.GENE_ID_COLNAME).str.contains(r"-[53]p$"))
        .then(pl.col(config.GENE_ID_COLNAME).str.extract(r"^(.*?)-[53]p$"))
        .otherwise(pl.col(config.GENE_ID_COLNAME))
        .alias(config.GENE_ID_COLNAME)
    )


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    logger.info(f"Converting IDs for count file {args.count_file.name}...")

    #############################################################
    # PARSING FILES
    #############################################################

    df = parse_count_table(args.count_file)

    if df.is_empty():
        msg = "COUNT FILE IS EMPTY"
        logger.warning(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)

    try:
        df = clean_ensembl_gene_id_versioning(df)
        df = clean_mirna_ids(df)
    except Exception as e:
        msg = f"ERROR CLEANING IDS in count file {args.count_file.name}: {e}"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)

    #############################################################
    # WRITING RESULTS
    #############################################################

    logger.info("Writing cleaned IDs")
    gene_ids_outfile = args.count_file.with_name(
        args.count_file.stem + CLEANED_GENE_IDS_SUFFIX
    )
    gene_ids = (
        df.select(config.GENE_ID_COLNAME)
        .sort(config.GENE_ID_COLNAME)
        .to_series()
        .to_list()
    )

    with open(gene_ids_outfile, "w") as fout:
        fout.write("\n".join(gene_ids))

    logger.info("Writing count file with cleaned IDs")
    count_outfile = args.count_file.with_name(
        args.count_file.stem + CLEANED_COUNTS_SUFFIX
    )
    df.write_parquet(count_outfile)


if __name__ == "__main__":
    main()
