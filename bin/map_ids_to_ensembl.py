#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import pandas as pd
from pathlib import Path
import argparse
import logging
import sys

from gprofiler_utils import convert_ids

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

RENAMED_FILE_SUFFIX = ".renamed.csv"
METADATA_FILE_SUFFIX = ".metadata.csv"
MAPPING_FILE_SUFFIX = ".mapping.csv"

ORIGINAL_GENE_ID_COLNAME = "original_gene_id"
ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"

##################################################################
# FUNCTIONS
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Map IDs to Ensembl")
    parser.add_argument(
        "--count-file", type=Path, required=True, help="Input file containing counts"
    )
    parser.add_argument(
        "--species", type=str, required=True, help="Species to convert IDs for"
    )
    parser.add_argument(
        "--custom-mappings", type=str, help="Optional file containing custom mappings"
    )
    return parser.parse_args()


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    count_file = args.count_file
    logger.info(
        f"Converting IDs for species {args.species} and count file {count_file.name}..."
    )

    #############################################################"
    # PARSING FILES
    #############################################################
    df = pd.read_csv(count_file, header=0, index_col=0)
    if df.empty:
        logger.error("Count file is empty! Aborting ID mapping...")
        sys.exit(100)

    df.index = df.index.astype(str)
    gene_ids = df.index.tolist()

    custom_mappings_dict = {}
    custom_mapping_file = args.custom_mappings
    if custom_mapping_file:
        if Path(custom_mapping_file).is_file():
            custom_mapping_df = pd.read_csv(custom_mapping_file)
            custom_mappings_dict = custom_mapping_df.set_index(
                ORIGINAL_GENE_ID_COLNAME
            )[ENSEMBL_GENE_ID_COLNAME].to_dict()

    gene_ids_left_to_map = [
        gene_id for gene_id in gene_ids
        if gene_id not in custom_mappings_dict.keys()
    ]
    logger.info(f"Number of genes left to map: {len(gene_ids_left_to_map)}")

    mapping_dict = {}
    gene_metadata_dfs = []

    #############################################################
    # QUERYING g:PROFILER SERVER
    #############################################################

    if gene_ids_left_to_map:
        mapping_dict, gene_metadata_dfs = convert_ids(gene_ids_left_to_map, args.species)

    # adding custom mappings
    mapping_dict.update(custom_mappings_dict)
    # if mapping dict is empty
    if not mapping_dict:
        logger.error(
            f"No mapping found for gene names in count file {count_file.name} "
            f"and for species {args.species}! "
            f"Example of gene names found in the provided dataframe: {df.index[:5].tolist()}"
            f"Count file is empty! Aborting ID mapping..."
        )
        sys.exit(101)

    #############################################################"
    # MAPPING GENE IDS IN DATAFRAME
    #############################################################
    # filtering the DataFrame to keep only the rows where the index can be mapped
    df = df.loc[df.index.isin(mapping_dict)]

    # renaming gene names to mapped ids using mapping dict
    df.index = df.index.map(mapping_dict)
    df.reset_index(inplace=True)
    df.rename(columns={"index": ENSEMBL_GENE_ID_COLNAME}, inplace=True)

    # TODO: check is there is another way to avoid duplicate gene names
    # sometimes different gene names have the same ensembl ID
    # for now, we just get the mean of values, but this is not ideal
    df = df.groupby(ENSEMBL_GENE_ID_COLNAME, as_index=False).mean()

    #############################################################"
    # WRITING OUTFILES
    #############################################################
    # writing to output file
    outfile = count_file.with_name(count_file.stem + RENAMED_FILE_SUFFIX)
    df.to_csv(outfile, index=False, header=True)

    # concatenating all metadata and ensuring there are no duplicates
    if gene_metadata_dfs:
        gene_metadata_df = pd.concat(gene_metadata_dfs, ignore_index=True)
        gene_metadata_df.drop_duplicates(inplace=True)
        # writing gene metadata to file
        metadata_file = count_file.with_name(count_file.stem + METADATA_FILE_SUFFIX)
        gene_metadata_df.to_csv(metadata_file, index=False, header=True)

    # making dataframe for mapping (only two columns: original and new)
    mapping_df = (
        pd.DataFrame(mapping_dict, index=[0])
        .T.reset_index()  # transpose: setting keys as indexes instead of columns
        .rename(columns={"index": ORIGINAL_GENE_ID_COLNAME, 0: ENSEMBL_GENE_ID_COLNAME})
    )
    mapping_file = count_file.with_name(count_file.stem + MAPPING_FILE_SUFFIX)
    mapping_df.to_csv(mapping_file, index=False, header=True)


if __name__ == "__main__":
    main()
