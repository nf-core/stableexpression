#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import json
import logging
from pathlib import Path

import config
import pandas as pd
import requests
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)
from tqdm.contrib.concurrent import process_map

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

GENE_IDS_CHUNKSIZE = 50  # max allowed by Ensembl REST API

ENSEMBL_REST_SERVER = "https://rest.ensembl.org"
SEQUENCE_INFO_EXT = "/sequence/id"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}
STOP_RETRY_AFTER_DELAY = 600

OUTFILE = "gene_ids_lengths.csv"


##################################################################
##################################################################
# FUNCTIONS
##################################################################
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Get GEO Datasets accessions")
    parser.add_argument(
        "--genes",
        type=Path,
        dest="gene_file",
        required=True,
        help="File containing gene IDs",
    )
    return parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUERIES TO ENSEMBL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_post_request_to_ensembl(gene_ids: list[str]) -> list[dict]:
    data = {"ids": gene_ids, "type": "cdna"}
    url = ENSEMBL_REST_SERVER + SEQUENCE_INFO_EXT
    response = requests.post(url, headers=HEADERS, data=json.dumps(data))
    if response.status_code == 200:
        response.raise_for_status()
    else:
        raise RuntimeError(
            f"Failed to retrieve data: encountered error {response.status_code}"
        )
    return response.json()


def get_gene_lengths(gene_ids: list[str]) -> list[dict]:
    records = send_post_request_to_ensembl(gene_ids)
    return [
        {
            config.GENE_ID_COLNAME: record["query"],
            config.CDNA_LENGTH_COLNAME: len(record["seq"]),
        }
        for record in records
        if record.get("query") is not None and record.get("seq") is not None
    ]


def chunk_list(lst: list, chunksize: int) -> list:
    """Splits a list into chunks of a given size.

    Args:
        lst (list): The list to split.
        chunksize (int): The size of each chunk.

    Returns:
        list: A list of chunks, where each chunk is a list of len(chunksize).
    """
    return [lst[i : i + chunksize] for i in range(0, len(lst), chunksize)]


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    args = parse_args()

    with open(args.gene_file, "r") as fin:
        gene_ids = [line.strip() for line in fin]

    gene_id_chunks = chunk_list(gene_ids, GENE_IDS_CHUNKSIZE)
    # getting gene lengths chunk by chunk
    records_list = process_map(get_gene_lengths, gene_id_chunks, max_workers=12)
    # flattening  list of lists into a single list
    records = [record for sublist in records_list for record in sublist]

    df = pd.DataFrame.from_dict(records)
    # taking the length of the longest transcript for each gene
    df = df.groupby(config.GENE_ID_COLNAME, as_index=False).agg(
        {config.CDNA_LENGTH_COLNAME: "max"}
    )

    df.to_csv(OUTFILE, index=False, header=True)


if __name__ == "__main__":
    main()
