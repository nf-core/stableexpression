#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

RENAMED_FILE_SUFFIX = ".renamed.csv"
METADATA_FILE_SUFFIX = ".metadata.csv"
MAPPING_FILE_SUFFIX = ".mapping.csv"

WARNING_REASON_FILE = "warning_reason.txt"
FAILURE_REASON_FILE = "failure_reason.txt"

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
        "--custom-mappings",
        type=Path,
        dest="custom_mapping_file",
        help="Optional file containing custom mappings",
    )
    return parser.parse_args()


def parse_table(file: Path, **kwargs):
    if file.suffix == ".csv":
        return pd.read_csv(file, header=0, **kwargs)
    else:  # .tsv
        return pd.read_csv(file, header=0, sep="\t", **kwargs)


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    logger.info(f"Converting IDs for count file {args.count_file.name}...")

    #############################################################
    # PARSING FILES
    #############################################################

    # whatever the name of the first col, rename it to "gene_id"
    df = parse_table(args.count_file, index_col=0)
    df.index.rename(config.GENE_ID_COLNAME, inplace=True)

    if df.empty:
        msg = "COUNT FILE IS EMPTY"
        logger.warning(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)

    df.index = df.index.astype(str)

    #############################################################
    # GETTING MAPPINGS
    #############################################################

    mapping_dict = {}
    if args.mapping_file is not None:
        mapping_df = parse_table(args.mapping_file)
        mapping_dict = mapping_df.set_index(config.ORIGINAL_GENE_ID_COLNAME)[
            config.GENE_ID_COLNAME
        ].to_dict()

    custom_mapping_dict = {}
    if args.custom_mapping_file is not None:
        custom_mapping_df = parse_table(args.custom_mapping_file)
        custom_mapping_dict = custom_mapping_df.set_index(
            config.ORIGINAL_GENE_ID_COLNAME
        )[config.GENE_ID_COLNAME].to_dict()

    mapping_dict |= custom_mapping_dict

    if not mapping_dict:
        raise ValueError("No mapping found")  # should not happen

    #############################################################
    # MAPPING GENE IDS IN DATAFRAME
    #############################################################

    # IMPORTANT: KEEPING ONLY GENES THAT HAVE BEEN CONVERTED
    # filtering the DataFrame to keep only the rows where the index can be mapped
    original_nb_genes = len(df)

    df = df.loc[df.index.isin(mapping_dict)]
    if df.empty:
        msg = "NO GENES WERE MAPPED"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)

    if len(df) < original_nb_genes:
        msg = f"Only {len(df) / original_nb_genes:.2%} of genes were mapped ({len(df)} out of {original_nb_genes})"
        logger.warning(msg)
        with open(WARNING_REASON_FILE, "a") as f:
            f.write(msg)
    else:
        logger.info(f"All genes were mapped ({len(df)} out of {original_nb_genes})")

    # renaming gene names to mapped ids using mapping dict
    df.index = df.index.map(mapping_dict)
    df.reset_index(inplace=True)

    # TODO: check is there is another way to avoid duplicate gene names
    # sometimes different gene names have the same Gene ID
    # for now, we just get the mean of values, but this is not ideal

    #############################################################
    # GENE COUNT HANDLING
    #############################################################

    # handling cases where multiple genes have the same Gene ID
    # since subsequent steps in the pipeline require integer values,
    # we need to ensure that the resulting DataFrame has integer values
    df = df.groupby(config.GENE_ID_COLNAME, as_index=False, sort=False).agg(
        lambda x: x.mean().astype(int)
    )

    #############################################################
    # WRITING OUTFILES
    #############################################################
    # writing to output file
    outfile = args.count_file.with_name(args.count_file.stem + RENAMED_FILE_SUFFIX)
    df.to_csv(outfile, index=False, header=True)

    # making dataframe for mapping (only two columns: original and new)
    mapping_df = (
        pd.DataFrame(mapping_dict, index=[0])
        .T.reset_index()  # transpose: setting keys as indexes instead of columns
        .rename(
            columns={
                "index": config.ORIGINAL_GENE_ID_COLNAME,
                0: config.GENE_ID_COLNAME,
            }
        )
    )
    mapping_file = args.count_file.with_name(args.count_file.stem + MAPPING_FILE_SUFFIX)
    mapping_df.to_csv(mapping_file, index=False, header=True)


if __name__ == "__main__":
    main()
