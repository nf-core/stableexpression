#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import pandas as pd
from pathlib import Path
import argparse
import logging
import sys

from gprofiler_utils import convert_ids
import config

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

RENAMED_FILE_SUFFIX = ".renamed.csv"
METADATA_FILE_SUFFIX = ".metadata.csv"
MAPPING_FILE_SUFFIX = ".mapping.csv"

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
        "--custom-mappings",
        type=str,
        dest="custom_mappings",
        help="Optional file containing custom mappings",
    )
    parser.add_argument(
        "--custom-metadata",
        type=str,
        dest="custom_metadata",
        help="Optional file containing custom metadata",
    )
    return parser.parse_args()


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    count_file = args.count_file
    custom_mapping_file = args.custom_mappings
    custom_metadata_file = args.custom_metadata

    logger.info(
        f"Converting IDs for species {args.species} and count file {count_file.name}..."
    )

    #############################################################"
    # PARSING FILES
    #############################################################
    df = pd.read_csv(count_file, header=0, index_col=0)
    if df.empty:
        logger.warning("Count file is empty! Aborting ID mapping...")
        sys.exit(0)

    df.index = df.index.astype(str)
    gene_ids = df.index.tolist()

    custom_mappings_dict = {}
    if custom_mapping_file:
        custom_mapping_df = pd.read_csv(custom_mapping_file)
        custom_mappings_dict = custom_mapping_df.set_index(
            config.ORIGINAL_GENE_ID_COLNAME
        )[config.ENSEMBL_GENE_ID_COLNAME].to_dict()

    gene_ids_left_to_map = [
        gene_id for gene_id in gene_ids if gene_id not in custom_mappings_dict
    ]
    logger.info(f"Number of genes left to map: {len(gene_ids_left_to_map)}")

    gene_metadata_dfs = []

    #############################################################
    # QUERYING g:PROFILER SERVER
    #############################################################

    if gene_ids_left_to_map:
        gprofiler_mapping_dict, gene_metadata_dfs = convert_ids(
            gene_ids_left_to_map, args.species
        )

    # overall mappings is the custom_mappings_dict complemented with gprofiler_mapping_dict
    mapping_dict = custom_mappings_dict | gprofiler_mapping_dict

    # if mapping dict is empty
    if not mapping_dict:
        logger.warning(
            f"No mapping found for gene names in count file {count_file.name} "
            f"and for species {args.species}! "
            f"Example of gene names found in the provided dataframe: {df.index[:5].tolist()}"
            f"Count file is empty! Aborting ID mapping..."
        )
        sys.exit(0)

    #############################################################"
    # MAPPING GENE IDS IN DATAFRAME
    #############################################################

    # IMPORTANT: KEEPING ONLY GENES THAT HAVE BEEN CONVERTED
    # filtering the DataFrame to keep only the rows where the index can be mapped
    df = df.loc[df.index.isin(mapping_dict)]

    # renaming gene names to mapped ids using mapping dict
    df.index = df.index.map(mapping_dict)
    df.reset_index(inplace=True)
    df.rename(columns={"index": config.ENSEMBL_GENE_ID_COLNAME}, inplace=True)

    # TODO: check is there is another way to avoid duplicate gene names
    # sometimes different gene names have the same ensembl ID
    # for now, we just get the mean of values, but this is not ideal
    df = df.groupby(config.ENSEMBL_GENE_ID_COLNAME, as_index=False).mean()

    #############################################################"
    # WRITING OUTFILES
    #############################################################
    # writing to output file
    outfile = count_file.with_name(count_file.stem + RENAMED_FILE_SUFFIX)
    df.to_csv(outfile, index=False, header=True)

    # if the user provides custom metadata file
    if custom_metadata_file:
        custom_metadata_df = pd.read_csv(custom_metadata_file)
        # prepending custom metadata in gene metadata
        gene_metadata_dfs = [custom_metadata_df] + gene_metadata_dfs

    # concatenating all metadata and ensuring there are no duplicates
    if gene_metadata_dfs:
        gene_metadata_df = pd.concat(gene_metadata_dfs, ignore_index=True)
        # dropping duplicates and keeping the first occurence
        gene_metadata_df.drop_duplicates(
            inplace=True, subset=[config.ENSEMBL_GENE_ID_COLNAME], keep="first"
        )
        # writing gene metadata to file
        metadata_file = count_file.with_name(count_file.stem + METADATA_FILE_SUFFIX)
        gene_metadata_df.to_csv(metadata_file, index=False, header=True)

    # making dataframe for mapping (only two columns: original and new)
    mapping_df = (
        pd.DataFrame(mapping_dict, index=[0])
        .T.reset_index()  # transpose: setting keys as indexes instead of columns
        .rename(
            columns={
                "index": config.ORIGINAL_GENE_ID_COLNAME,
                0: config.ENSEMBL_GENE_ID_COLNAME,
            }
        )
    )
    mapping_file = count_file.with_name(count_file.stem + MAPPING_FILE_SUFFIX)
    mapping_df.to_csv(mapping_file, index=False, header=True)


if __name__ == "__main__":
    main()
