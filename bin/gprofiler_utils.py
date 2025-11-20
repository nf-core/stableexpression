#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging

import config
import pandas as pd
import requests
from requests.exceptions import ConnectionError, HTTPError
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

GPROFILER_CONVERT_API_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
GPROFILER_CONVERT_BETA_API_ENDPOINT = (
    "https://biit.cs.ut.ee/gprofiler_beta/api/convert/convert/"
)

CHUNKSIZE = 2000  # number of IDs to convert at a time - may create trouble if > 2000

COLS_TO_KEEP = ["incoming", "converted", "name", "description"]
DESCRIPTION_PART_TO_REMOVE_REGEX = r"\s*\[Source:.*?\]"

GPROFILER_ERROR_MESSAGE = (
    "g:Profiler servers (main and beta) seem to be down... Please retry later... "
    "If you have gene ID mappings and / or gene metadata for these datasets, you can provide them "
    "directly using the `--gene_id_mapping` and `--gene_metadata` parameters respectively, "
    "and by skipping the g:Profiler ID mapping step with `--skip_id_mapping`."
)


##################################################################
# FUNCTIONS
##################################################################


class GProfilerConnectionError(Exception):
    pass


def format_species_name(species: str):
    """
    Format a species name into a format accepted by g:Profiler.
    Example: Arabidopsis thaliana -> athaliana

    Parameters
    ----------
    species : str
        The species name.

    Returns
    -------
    str
        The formatted species name.
    """
    splitted_species = species.lower().replace("_", " ").split(" ")
    return splitted_species[0][0] + splitted_species[1]


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def request_conversion(
    gene_ids: list,
    species: str,
    target_database: str,
    url: str = GPROFILER_CONVERT_API_ENDPOINT,
    attempts: int = 0,
) -> list[str]:
    """
    Send a request to the g:Profiler API to convert a list of gene IDs.

    Parameters
    ----------
    gene_ids : list
        The list of gene IDs to convert.
    species : str
        The species to convert the IDs for.
    url : str, optionalrequest_conversion
        The URL to send the request to, by default GPROFILER_CONVERT_API_ENDPOINT
    attempts : int, optional
        The number of attempts already performed, by default 0

    Returns
    -------
    list
        The list of dicts corresponding to the converted IDs.
    """

    # formatting species for g:Profiler
    organism = format_species_name(species)

    if attempts > 0:
        logger.warning(
            "g:Profiler main server appears down, trying with the beta server..."
        )

    server_appears_down = False

    try:
        response = requests.post(
            url=url,
            json={"organism": organism, "query": gene_ids, "target": target_database},
        )
    except requests.exceptions.ConnectionError:
        server_appears_down = True
    else:
        try:
            response.raise_for_status()
        except (HTTPError, ConnectionError) as err:
            if str(response.status_code).startswith("5"):  # error 500 -> 509
                server_appears_down = True
            else:
                logger.error(
                    f"Error {response.status_code} while converting IDs: {err}"
                )
                raise err

    if server_appears_down:
        if attempts == 0:
            logger.warning(
                "g:Profiler main server appears down, trying with the beta server..."
            )
            return request_conversion(
                gene_ids,
                species,
                target_database=target_database,
                url=GPROFILER_CONVERT_BETA_API_ENDPOINT,  # backup endpoint
                attempts=1,
            )
        else:
            # both servers appear down, we stop here...
            logger.error(GPROFILER_ERROR_MESSAGE)
            raise GProfilerConnectionError(GPROFILER_ERROR_MESSAGE)

    else:
        return response.json()["result"]


def convert_chunk_of_ids(
    gene_ids: list, species: str, gprofiler_target_db: str
) -> tuple[dict, pd.DataFrame]:
    """
    Wrapper function that converts a list of gene IDs to another namespace.

    Parameters
    ----------
    species : str
        The species to convert the IDs for.
    gene_ids : list
        The IDs to convert.
    target_database : str
        The target database to convert to.

    Returns
    -------
    dict
        A dictionary where the keys are the original IDs and the values are the converted IDs.
    """

    results = request_conversion(gene_ids, species, gprofiler_target_db)
    df = pd.DataFrame.from_records(results)

    if df.empty:
        return {}, pd.DataFrame()

    # keeping only rows where 'converted' is not null and only the columns of interest
    df = df.loc[df["converted"] != "None", COLS_TO_KEEP]

    # dict associating incoming IDs to converted IDs
    mapping_dict = df.set_index("incoming").to_dict()["converted"]

    # DataFrame associating converted IDs to name and description
    meta_df = df.drop(columns=["incoming"]).rename(
        columns={"converted": config.GENE_ID_COLNAME}
    )

    meta_df["name"] = meta_df["name"].str.replace(",", ";")

    # Extract the part before '[Source:...]', or the whole string if not found
    meta_df["description"] = (
        meta_df["description"]
        .str.replace(DESCRIPTION_PART_TO_REMOVE_REGEX, "", regex=True)
        .str.replace(",", ";")
    )

    return mapping_dict, meta_df


def chunk_list(lst: list, chunksize: int) -> list:
    """Splits a list into chunks of a given size.

    Args:
        lst (list): The list to split.
        chunksize (int): The size of each chunk.

    Returns:
        list: A list of chunks, where each chunk is a list of len(chunksize).
    """
    return [lst[i : i + chunksize] for i in range(0, len(lst), chunksize)]


def convert_ids(
    ids: list[str], species: str, gprofiler_target_db: str
) -> tuple[dict, pd.DataFrame]:
    mapping_dict = {}
    gene_metadata_dfs = []

    chunks = chunk_list(ids, chunksize=CHUNKSIZE)
    for chunk_gene_ids in chunks:
        # converting to Gene IDs for all IDs comprised in this chunk
        gene_mapping, meta_df = convert_chunk_of_ids(
            chunk_gene_ids, species, gprofiler_target_db
        )
        mapping_dict.update(gene_mapping)
        gene_metadata_dfs.append(meta_df)

    return mapping_dict, gene_metadata_dfs
