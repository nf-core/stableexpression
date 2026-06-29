#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
from collections import Counter

import httpx
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ALL_EXP_URL = "https://www.ebi.ac.uk/gxa/json/experiments/"

##################################################################
##################################################################
# FUNCTIONS
##################################################################
##################################################################


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def get_data(url: str) -> dict:
    """
    Queries a URL and returns the data as a JSON object

    Parameters
    ----------
    url : str
        The URL to query

    Returns
    -------
    data : dict
        The JSON object returned by the query

    Raises
    ------
    RuntimeError
        If the query fails
    """
    response = httpx.get(url)
    response.raise_for_status()
    return response.json()


def get_eatlas_experiments():
    """
    Gets all experiments from Expression Atlas

    Parameters
    ----------

    Returns
    -------
    experiments : list
        A list of experiment dictionaries
    """
    data = get_data(ALL_EXP_URL)
    return data["experiments"]


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    experiments = get_eatlas_experiments()
    species = Counter(
        sorted([" ".join(exp["species"].split(" ")[:2]) for exp in experiments])
    )
    print(species)
    print(len(species))


if __name__ == "__main__":
    main()
