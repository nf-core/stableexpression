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

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

GENE_IDS_CHUNKSIZE = 50  # max allowed by Ensembl REST API

ENSEMBL_REST_SERVER = "https://rest.ensembl.org/"
SPECIES_INFO_EXT = "info/genomes/taxonomy/{species}"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}
STOP_RETRY_AFTER_DELAY = 600


##################################################################
##################################################################
# FUNCTIONS
##################################################################
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Get GEO Datasets accessions")
    parser.add_argument(
        "--species",
        type=str,
        dest="species",
        required=True,
        help="Species name",
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
def send_get_request_to_ensembl(url: str) -> list[dict]:
    logger.info(f"Sending GET request to {url}")
    response = requests.get(url, headers=HEADERS)
    if response.status_code == 200:
        response.raise_for_status()
    else:
        raise RuntimeError(
            f"Failed to retrieve data: encountered error {response.status_code}"
        )
    return response.json()


def get_species_division(species: str) -> str:
    url = ENSEMBL_REST_SERVER + SPECIES_INFO_EXT.format(species=species)
    data = send_get_request_to_ensembl(url)
    if len(data) == 0:
        raise ValueError(f"No division found for {species}")
    elif len(data) > 1:
        logger.warning(
            f"Multiple divisions found for {species}. Keeping the first one."
        )
    return data[0]["division"]


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    args = parse_args()

    species_division = get_species_division(args.species)
    print(species_division)


if __name__ == "__main__":
    main()
