#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import pandas as pd
from gprofiler_utils import convert_ids

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
    parser.add_argument(
        "--memory", type=str, dest="memory", required=True, help="Memory in GB"
    )
    return parser.parse_args()


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    with open(args.gene_id_file, "r") as fin:
        gene_ids = list(set([line.strip() for line in fin]))

    logger.info(f"Converting {len(gene_ids)} IDs for species {args.species} ")

    #############################################################
    # QUERYING g:PROFILER SERVER
    #############################################################

    gene_metadata_dfs = []

    mapping_dict, gene_metadata_dfs = convert_ids(
        gene_ids, args.species, args.gprofiler_target_db
    )

    if not mapping_dict:
        msg = (
            f"No mapping found for gene IDs such as {' '.join(gene_ids[:5])} on species {args.species} "
            + f"and g:Profiler target database {args.gprofiler_target_db}"
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

    gene_metadata_df = pd.concat(gene_metadata_dfs, ignore_index=True)
    # dropping duplicates and keeping the first occurence
    gene_metadata_df.drop_duplicates(
        subset=[config.GENE_ID_COLNAME], keep="first"
    ).sort_values(by=config.GENE_ID_COLNAME).to_csv(
        METADATA_OUTFILE, index=False, header=True
    )


if __name__ == "__main__":
    main()
