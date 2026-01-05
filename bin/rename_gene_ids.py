#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import polars as pl
from common import parse_count_table, parse_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

RENAMED_FILE_SUFFIX = ".renamed.parquet"

WARNING_REASON_FILE = "warning_reason.txt"
FAILURE_REASON_FILE = "failure_reason.txt"

UNMAPPED_FILE_SUFFIX = "unmapped.txt"
NOT_VALID_FILE_SUFFIX = "not_valid.txt"
MERGED_FILE_SUFFIX = "merged.txt"
FINAL_FILE_SUFFIX = "final.txt"

##################################################################
# FUNCTIONS
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Rename gene IDs using mapped IDs")
    parser.add_argument(
        "--count-file", type=Path, required=True, help="Input file containing counts"
    )
    parser.add_argument(
        "--mappings",
        type=Path,
        dest="mapping_file",
        help="Mapping file containing gene IDs",
    )
    parser.add_argument(
        "--valid-gene-ids",
        type=Path,
        dest="valid_gene_ids_file",
        help="File containing valid gene IDs",
    )
    return parser.parse_args()


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

    #############################################################
    # GETTING MAPPINGS
    #############################################################

    mapping_df = parse_table(args.mapping_file)
    mapping_dict = dict(
        zip(
            mapping_df[config.ORIGINAL_GENE_ID_COLNAME],
            mapping_df[config.GENE_ID_COLNAME],
        )
    )

    #############################################################
    # MAPPING GENE IDS IN DATAFRAME
    #############################################################

    # IMPORTANT: KEEPING ONLY GENES THAT HAVE BEEN CONVERTED
    # filtering the DataFrame to keep only the rows where the index can be mapped
    original_nb_genes = len(df)

    rejected_df = df.filter(~pl.col(config.GENE_ID_COLNAME).is_in(mapping_dict.keys()))
    nb_unmapped_genes = len(rejected_df)

    # df = df.loc[df.index.isin(mapping_dict)]
    df = df.filter(pl.col(config.GENE_ID_COLNAME).is_in(mapping_dict.keys()))
    nb_mapped_genes = len(df)

    with open(UNMAPPED_FILE_SUFFIX, "w") as f:
        f.write(str(nb_unmapped_genes))

    if df.is_empty():
        sample_size = min(5, nb_unmapped_genes)
        example_rejected_genes = (
            rejected_df[config.GENE_ID_COLNAME].head(sample_size).to_list()
        )
        msg = f"NO GENES WERE MAPPED. EXAMPLE OF GENE IDS: {example_rejected_genes}"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)

        with open(NOT_VALID_FILE_SUFFIX, "w") as f:
            f.write("0")
        with open(MERGED_FILE_SUFFIX, "w") as f:
            f.write("0")
        with open(FINAL_FILE_SUFFIX, "w") as f:
            f.write("0")

        sys.exit(0)

    if len(df) < original_nb_genes:
        sample_size = min(5, nb_unmapped_genes)
        example_rejected_genes = (
            rejected_df[config.GENE_ID_COLNAME].head(sample_size).to_list()
        )
        msg = (
            f"{nb_mapped_genes / original_nb_genes:.2%} of genes were mapped ({nb_mapped_genes} out of {original_nb_genes}). "
            + f"Example of unmapped genes: {example_rejected_genes}"
        )
        logger.warning(msg)
        with open(WARNING_REASON_FILE, "a") as f:
            f.write(msg)
    else:
        logger.info(
            f"All genes were mapped ({nb_mapped_genes} out of {original_nb_genes})"
        )

    logger.info("Renaming gene names")
    # renaming gene names to mapped ids using mapping dict
    df = df.with_columns(
        pl.col(config.GENE_ID_COLNAME)
        .replace(mapping_dict)
        .alias(config.GENE_ID_COLNAME)
    )

    #############################################################
    # GETTING VALID GENE IDS
    #############################################################

    logger.info("Keeping only genes with sufficient occurrence over datasets")
    nb_genes_before_validation = len(df)

    with open(args.valid_gene_ids_file, "r") as fin:
        valid_gene_ids = [line.strip() for line in fin.readlines()]

    df = df.filter(pl.col(config.GENE_ID_COLNAME).is_in(valid_gene_ids))

    nb_not_valid_genes = nb_genes_before_validation - len(df)
    logger.info(
        f"{nb_not_valid_genes} ({nb_not_valid_genes / nb_genes_before_validation:.2%}) genes were not valid"
    )

    with open(NOT_VALID_FILE_SUFFIX, "w") as f:
        f.write(str(nb_not_valid_genes))

    if df.is_empty():
        msg = "NO GENES LEFT AFTER REMOVING RARE GENE IDS"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)

        with open(MERGED_FILE_SUFFIX, "w") as f:
            f.write("0")
        with open(FINAL_FILE_SUFFIX, "w") as f:
            f.write("0")

        sys.exit(0)

    #############################################################
    # GENE COUNT HANDLING
    #############################################################

    # handling cases where multiple genes have the same Gene ID
    # since subsequent steps in the pipeline require integer values,
    # we need to ensure that the resulting DataFrame has integer values

    # TODO: check is there is another way to avoid duplicate gene names
    # sometimes different gene names have the same Gene ID
    # for now, we just get the mean of values, but this is not ideal

    logger.info("Computing mean counts for genes with duplicate IDs")
    df = df.group_by(config.GENE_ID_COLNAME, maintain_order=True).agg(
        pl.exclude(config.GENE_ID_COLNAME).mean()
    )

    #############################################################
    # WRITING OUTFILES
    #############################################################

    nb_merged = nb_mapped_genes - len(df)
    with open(MERGED_FILE_SUFFIX, "w") as f:
        f.write(str(nb_merged))
    with open(FINAL_FILE_SUFFIX, "w") as f:
        f.write(str(len(df)))

    logger.info("Writing output file")
    outfilename = args.count_file.with_suffix(RENAMED_FILE_SUFFIX).name
    df.write_parquet(outfilename)


if __name__ == "__main__":
    main()
