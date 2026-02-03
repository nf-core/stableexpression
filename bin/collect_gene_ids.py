#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from collections import Counter
from pathlib import Path

from tqdm import tqdm

import config

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

UNIQUE_GENE_IDS_OUTFILE = "unique_gene_ids.txt"
GENE_ID_OCCURRENCES_OUTFILE = "gene_id_occurrences.csv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Collect gene IDs from count files")
    parser.add_argument(
        "--ids", type=str, dest="gene_id_files", required=True, help="Gene ID files"
    )
    parser.add_argument(
        "--cpus", type=int, dest="nb_cpus", required=True, help="Number of CPUs"
    )
    parser.add_argument(
        "--memory", type=str, dest="memory", required=True, help="Memory in GB"
    )
    return parser.parse_args()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    gene_id_files = [Path(file) for file in args.gene_id_files.split(" ")]
    logger.info(f"Getting gene IDs from {len(gene_id_files)} files")

    unique_gene_ids = set()
    counter = Counter()
    for gene_id_file in tqdm(gene_id_files):
        with open(gene_id_file, "r") as fin:
            gene_ids = [line.strip() for line in fin]
            unique_gene_ids.update(gene_ids)
            counter.update(gene_ids)

    with open(UNIQUE_GENE_IDS_OUTFILE, "w") as fout:
        fout.write("\n".join([str(gene_id) for gene_id in sorted(unique_gene_ids)]))

    with open(GENE_ID_OCCURRENCES_OUTFILE, "w") as fout:
        fout.write(
            f"{config.ORIGINAL_GENE_ID_COLNAME},{config.GENE_ID_COUNT_COLNAME}\n"
        )
        for gene_id, count in sorted(counter.items()):
            fout.write(f"{gene_id},{count}\n")


if __name__ == "__main__":
    main()
