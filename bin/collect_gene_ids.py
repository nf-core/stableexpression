#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import pandas as pd
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ALL_GENE_IDS_OUTFILE = "all_gene_ids.txt"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Collect gene IDs from count files")
    parser.add_argument(
        "--counts", type=str, dest="count_files", required=True, help="Count files"
    )
    return parser.parse_args()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def parse_table(file: Path):
    if file.suffix == ".csv":
        return pd.read_csv(file, header=0, index_col=0)
    else:  # .tsv
        return pd.read_csv(file, header=0, index_col=0, sep="\t")


def main():
    args = parse_args()
    count_files = [Path(file) for file in args.count_files.split(" ")]
    logger.info(f"Getting gene IDs from {len(count_files)} count files")

    all_gene_ids = set()
    for count_file in tqdm(count_files):
        df = parse_table(count_file)
        all_gene_ids.update(list(df.index))

    with open(ALL_GENE_IDS_OUTFILE, "w") as f:
        f.write("\n".join([str(gene_id) for gene_id in all_gene_ids]))


if __name__ == "__main__":
    main()
