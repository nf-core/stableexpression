#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import parse_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE = "valid_gene_ids.txt"

##################################################################
# FUNCTIONS
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Get genes with good occurrence")
    parser.add_argument(
        "--occurrences",
        type=Path,
        required=True,
        dest="gene_id_occurrence_file",
        help="Input file containing gene ID occurrences",
    )
    parser.add_argument(
        "--mappings",
        type=Path,
        required=True,
        dest="mapping_file",
        help="Mapping file containing gene IDs",
    )
    parser.add_argument(
        "--nb-datasets",
        type=int,
        required=True,
        dest="nb_datasets",
        help="Number of datasets",
    )
    parser.add_argument(
        "--min-freq-occurrence",
        type=float,
        required=True,
        dest="min_freq_occurrence",
        help="Minimum frequency of occurrences for a gene among all datasets",
    )
    return parser.parse_args()


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    # taking lower bound of threshold
    occurrence_threshold = int(args.nb_datasets * args.min_freq_occurrence)
    logger.info(
        f"Occurrence threshold: at least {occurrence_threshold} occurrence(s) among {args.nb_datasets} dataset(s)"
    )

    original_gene_id_occurrence_df = parse_table(args.gene_id_occurrence_file)
    mapping_df = parse_table(args.mapping_file)
    nb_mapped_genes = len(mapping_df)

    df = original_gene_id_occurrence_df.join(
        mapping_df,
        on=config.ORIGINAL_GENE_ID_COLNAME,
    )

    total_gene_id_occurrence_df = df.group_by(config.GENE_ID_COLNAME).agg(
        pl.col(config.GENE_ID_COUNT_COLNAME).sum().alias("total")
    )

    df = df.join(
        total_gene_id_occurrence_df,
        on=config.GENE_ID_COLNAME,
    ).filter(pl.col("total") >= occurrence_threshold)

    valid_gene_ids = df.select(config.GENE_ID_COLNAME).unique().to_series().to_list()

    with open(OUTFILE, "w") as f:
        f.write("\n".join(valid_gene_ids))

    nb_valid_genes = len(valid_gene_ids)

    logger.info(
        f"Found {nb_valid_genes} valid gene IDs ({nb_valid_genes / nb_mapped_genes:.2%})"
    )


if __name__ == "__main__":
    main()
