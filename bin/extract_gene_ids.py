#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import parse_count_table
from resource_management import set_max_memory

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

CLEANED_GENE_IDS_SUFFIX = ".gene_ids.txt"


def parse_args():
    parser = argparse.ArgumentParser("Rename gene IDs using mapped IDs")
    parser.add_argument(
        "--count-file", type=Path, required=True, help="Input file containing counts"
    )
    parser.add_argument(
        "--memory", type=str, dest="memory", required=True, help="Memory in GB"
    )
    return parser.parse_args()


def get_sorted_gene_ids(df: pl.DataFrame):
    return (
        df.select(config.GENE_ID_COLNAME)
        .sort(config.GENE_ID_COLNAME)
        .to_series()
        .to_list()
    )


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    set_max_memory(args.memory)

    logger.info(f"Converting IDs for count file {args.count_file.name}...")

    df = parse_count_table(args.count_file)

    logger.info("Writing cleaned IDs")
    gene_ids_outfile = args.count_file.with_name(
        args.count_file.stem + CLEANED_GENE_IDS_SUFFIX
    )
    gene_ids = get_sorted_gene_ids(df)

    with open(gene_ids_outfile, "w") as fout:
        fout.write("\n".join(gene_ids))


if __name__ == "__main__":
    main()
