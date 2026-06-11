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

VALID_GENE_IDS_OUTFILE = "valid_gene_ids.txt"
TOTAL_OCCURRENCES_OUTFILE = "total_gene_id_occurrence_quantiles.csv"

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
        "--min-occurrence-frequency",
        type=float,
        required=True,
        dest="min_occurrence_frequency",
        help="Minimum frequency of occurrences for a gene among all datasets",
    )
    parser.add_argument(
        "--min-occurrence-quantile",
        type=float,
        required=True,
        dest="min_occurrence_quantile",
        help="Minimum frequency of occurrences for a gene among all datasets",
    )
    return parser.parse_args()


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    original_gene_id_occurrence_df = parse_table(args.gene_id_occurrence_file)
    mapping_df = parse_table(args.mapping_file)
    nb_mapped_genes = len(mapping_df)

    df = original_gene_id_occurrence_df.join(
        mapping_df,
        on=config.ORIGINAL_GENE_ID_COLNAME,
    )

    total_gene_id_occurrence_df = df.group_by(config.GENE_ID_COLNAME).agg(
        pl.col(config.GENE_ID_COUNT_COLNAME).sum().alias("total_occurrences")
    )

    df = (
        df.join(
            total_gene_id_occurrence_df,
            on=config.GENE_ID_COLNAME,
        )
        .with_columns(
            total_occurrences_quantile=(
                pl.col("total_occurrences").rank(method="max")
                / pl.col("total_occurrences").count()
            ),
            total_occurrences_frequency=(
                pl.col("total_occurrences") / args.nb_datasets
            ),
        )
        .select(
            [
                config.GENE_ID_COLNAME,
                "total_occurrences_frequency",
                "total_occurrences_quantile",
            ]
        )
        .unique()
    )

    # sorting (for output consistency)
    df = df.sort(["total_occurrences_quantile", "gene_id"], descending=[True, False])

    # writing total occurrences in a csv before filtering
    df.select([config.GENE_ID_COLNAME, "total_occurrences_quantile"]).write_csv(
        TOTAL_OCCURRENCES_OUTFILE
    )

    # filtering genes
    valid_gene_ids = (
        df.filter(pl.col("total_occurrences_quantile") >= args.min_occurrence_quantile)
        .filter(pl.col("total_occurrences_frequency") >= args.min_occurrence_frequency)
        .select(config.GENE_ID_COLNAME)
        .unique()
        .to_series()
        .to_list()
    )

    with open(VALID_GENE_IDS_OUTFILE, "w") as f:
        f.write("\n".join(valid_gene_ids))

    nb_valid_genes = len(valid_gene_ids)

    logger.info(
        f"Found {nb_valid_genes} valid gene IDs ({nb_valid_genes / nb_mapped_genes:.2%})"
    )


if __name__ == "__main__":
    main()
