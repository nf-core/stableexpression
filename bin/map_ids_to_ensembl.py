#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import requests
import pandas as pd
from pathlib import Path
import argparse
import logging
import sys

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

RENAMED_FILE_SUFFIX = "_renamed.csv"
METADATA_FILE_SUFFIX = "_metadata.csv"
MAPPING_FILE_SUFFIX = "_mapping.csv"

CHUNKSIZE = 2000  # number of IDs to convert at a time - may create trouble if > 2000

GPROFILER_CONVERT_API_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
GPROFILER_CONVERT_BETA_API_ENDPOINT = (
    "https://biit.cs.ut.ee/gprofiler_beta/api/convert/convert/"
)

TARGET_DATABASE = "ENSG"  # Ensembl database
COLS_TO_KEEP = ["incoming", "converted", "name", "description"]
DESCRIPTION_PART_TO_REMOVE_REGEX = r"\s*\[Source:.*?\]"
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


def chunk_list(lst: list, chunksize: int):
    """Splits a list into chunks of a given size.

    Args:
        lst (list): The list to split.
        chunksize (int): The size of each chunk.

    Returns:
        list: A list of chunks, where each chunk is a list of len(chunksize).
    """
    return [lst[i : i + chunksize] for i in range(0, len(lst), chunksize)]


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
    url : str, optional
        The URL to send the request to, by default GPROFILER_CONVERT_API_ENDPOINT
    attempts : int, optional
        The number of attempts already performed, by default 0

    Returns
    -------
    list
        The list of dicts corresponding to the converted IDs.
    """
    response = requests.post(
        url=url,
        json={"organism": species, "query": gene_ids, "target": target_database},
    )
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as err:
        if err.response.status_code == 502:
            # server appears down
            if attempts == 0:
                # we only tried with the main server, we try with the beta server
                return request_conversion(
                    gene_ids,
                    species,
                    target_database=target_database,
                    url=GPROFILER_CONVERT_BETA_API_ENDPOINT,
                    attempts=1,
                )
            else:
                # both servers appear down, we stop here...
                logger.error(
                    "g:Profiler servers (main and beta) seem to be down... Please retry later... "
                    "If you have gene ID mappings and / or gene metadata for these datasets, you can provide them "
                    "directly using the `--gene_id_mapping` and `--gene_metadata` parameters respectively, "
                    "and by skipping the g:Profiler ID mapping step with `--skip_gprofiler`."
                )
                sys.exit(102)

        logger.error(f"Error {err.response.status_code} while converting IDs: {err}")
        sys.exit(101)

    return response.json()["result"]


def convert_ids(gene_ids: list, species: str):
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

    results = request_conversion(gene_ids, species, TARGET_DATABASE)
    df = pd.DataFrame.from_records(results)

    if df.empty:
        return {}

    # keeping only rows where 'converted' is not null and only the columns of interest
    df = df.loc[df["converted"] != "None", COLS_TO_KEEP]

    # dict associating incoming IDs to converted IDs
    mapping_dict = df.set_index("incoming").to_dict()["converted"]

    # DataFrame associating converted IDs to name and description
    meta_df = df.drop(columns=["incoming"]).rename(
        columns={"converted": ENSEMBL_GENE_ID_COLNAME}
    )

    meta_df["name"] = meta_df["name"].str.replace(",", ";")

    # Extract the part before '[Source:...]', or the whole string if not found
    meta_df["description"] = (
        meta_df["description"]
        .str.replace(DESCRIPTION_PART_TO_REMOVE_REGEX, "", regex=True)
        .str.replace(",", ";")
    )

    return mapping_dict, meta_df


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    count_file = args.count_file
    species_name = format_species_name(args.species)
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
        gene_id for gene_id in gene_ids if gene_id not in custom_mappings_dict.keys()
    ]
    logger.info(f"Number of genes left to map: {len(gene_ids_left_to_map)}")

    #############################################################
    # QUERYING g:PROFILER SERVER
    #############################################################
    mapping_dict = {}
    gene_metadata_dfs = []

    if gene_ids_left_to_map:
        chunks = chunk_list(gene_ids_left_to_map, chunksize=CHUNKSIZE)
        for chunk_gene_ids in chunks:
            # converting to Ensembl IDs for all IDs comprised in this chunk
            gene_mapping, meta_df = convert_ids(chunk_gene_ids, species_name)
            mapping_dict.update(gene_mapping)
            gene_metadata_dfs.append(meta_df)

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
