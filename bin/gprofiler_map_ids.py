#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import pandas as pd
from gprofiler_utils import convert_ids, get_candidate_organism_identifiers, format_species_name

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

MAPPED_GENE_IDS_OUTFILE = "mapped_gene_ids.csv"
METADATA_OUTFILE = "gene_metadata.csv"

TARGET_DATABASE_CHOICES = ["ENTREZGENE", "ENSG"]

FAILURE_REASON_FILE = "failure_reason.txt"

##################################################################
# FUNCTIONS
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Map IDs using g:Profiler")
    parser.add_argument(
        "--gene-ids",
        type=Path,
        dest="gene_id_file",
        required=True,
        help="Input file containing gene IDs",
    )
    parser.add_argument(
        "--species", type=str, required=True, help="Species to convert IDs for"
    )
    parser.add_argument(
        "--target-db",
        type=str,
        dest="gprofiler_target_db",
        required=True,
        choices=TARGET_DATABASE_CHOICES,
        help="Target database to convert IDs to",
    )
    return parser.parse_args()


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    with open(args.gene_id_file, "r") as fin:
        gene_ids = list({line.strip() for line in fin})

    logger.info(f"Converting {len(gene_ids)} IDs for species {args.species} ")

    #############################################################
    # QUERYING g:PROFILER SERVER FOR ALL POSSIBLE CANDIDATE ORGANISM IDENTIFIERS
    #############################################################

    candidate_organism_identifiers = get_candidate_organism_identifiers(args.species)

    if not candidate_organism_identifiers:
        raise ValueError(f"Species '{args.species}' not found in g:Profiler database.")

    #############################################################
    # CONVERT IDS BY CHUNKS FOR ALL CANDIDATE IDENTIFEIRS
    #############################################################

    candidate_mapping_dicts = []
    candidate_gene_metadata_df = []

    for organism in candidate_organism_identifiers:

        mapping_dict, gene_metadata_df = convert_ids(
            gene_ids, organism, args.gprofiler_target_db
        )
        candidate_mapping_dicts.append(mapping_dict)
        candidate_gene_metadata_df.append(gene_metadata_df)

    #############################################################
    # SELECTING THE BEST GENE ID MAPPING
    #############################################################

    # if multiple candidate identifiers
    # selecting the one that provides the best mapping with our gene IDs
    mapping_dict_lengths = [len(mapping_dict) for mapping_dict in candidate_mapping_dicts]
    # logging the depth of mapping for each candidate identifier
    for organism, mapping_dict_length in zip(candidate_organism_identifiers, mapping_dict_lengths):
        logger.info(f"Mapping size for {organism}: {mapping_dict_length}")

    best_mapping_index = mapping_dict_lengths.index(max(mapping_dict_lengths))
    mapping_dict = candidate_mapping_dicts[best_mapping_index]
    gene_metadata_df = candidate_gene_metadata_df[best_mapping_index]

    logger.info(f"Chosen organism: {candidate_organism_identifiers[best_mapping_index]}")

    if not mapping_dict: # no mapping of gene IDs, whatever the candidate identifier
        msg = (
            f"No mapping found for {args.species} against g:Profiler target database {args.gprofiler_target_db}. "
            + f"Example of unmapped gene IDs: {' '.join(gene_ids[:5])} "
        )
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as fout:
            fout.write(msg)
        sys.exit(100)

    #############################################################
    # WRITING MAPPING
    #############################################################

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
        .sort_values(by=config.ORIGINAL_GENE_ID_COLNAME)
    )
    mapping_df.to_csv(MAPPED_GENE_IDS_OUTFILE, index=False, header=True)

    #############################################################
    # WRITING METADATA
    #############################################################

    # dropping duplicates and keeping the first occurence
    gene_metadata_df.drop_duplicates(
        subset=[config.GENE_ID_COLNAME], keep="first"
    ).sort_values(by=config.GENE_ID_COLNAME).to_csv(
        METADATA_OUTFILE, index=False, header=True
    )


if __name__ == "__main__":
    main()
