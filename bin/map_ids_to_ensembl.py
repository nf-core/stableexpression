#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import requests
import pandas as pd
from pathlib import Path
import argparse
import json


class NoIDFoundException(Exception):
    pass


##################################################################
# CONSTANTS
##################################################################

RENAMED_FILE_SUFFIX = "_renamed.csv"
MAPPING_FILE_SUFFIX = "_mapping.json"
CHUNKSIZE = 2000  # number of IDs to convert at a time - may create trouble if > 2000

GPROFILER_CONVERT_API_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
TARGET_DATABASE = "ENSG"  # Ensembl database


##################################################################
# FUNCTIONS
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Map IDs to Ensembl")
    parser.add_argument("--count-file", type=Path, help="Input file containing counts")
    parser.add_argument("--species", type=str, help="Species to convert IDs for")
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
    gene_ids: list, species: str, target_database: str
) -> list[dict]:
    """
    Send a request to the g:Profiler API to convert a list of gene IDs.

    Parameters
    ----------
    gene_ids : list
        The list of gene IDs to convert.
    species : str
        The species to convert the IDs for.
    target_database : str
        The target database to convert to.

    Returns
    -------
    list
        The list of dicts corresponding to the converted IDs.
    """
    response = requests.post(
        url=GPROFILER_CONVERT_API_ENDPOINT,
        json={"organism": species, "query": gene_ids, "target": TARGET_DATABASE},
    )
    response.raise_for_status()
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
    df = df.loc[df["converted"] != "None", ["incoming", "converted"]]
    # changing index
    df.set_index("incoming", inplace=True)

    return df.to_dict()["converted"]


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    count_file = args.count_file
    species_name = format_species_name(args.species)
    print(
        f"Converting IDs for species {species_name} and count file {count_file.name}..."
    )

    df = pd.read_csv(count_file, header=0, index_col=0)
    df.index = df.index.astype(str)

    gene_ids = df.index.tolist()
    mapping_dict = {}

    chunks = chunk_list(gene_ids, chunksize=CHUNKSIZE)
    for chunk_gene_ids in chunks:
        # converting to Ensembl IDs for all IDs comprised in this chunk
        gene_mapping = convert_ids(chunk_gene_ids, species_name)
        mapping_dict.update(gene_mapping)

    if not mapping_dict:  # if mapping dict is empty
        raise NoIDFoundException(
            f"No mapping found for gene names in count file {count_file.name} "
            f"and for species {species_name}! "
            f"Example of gene names found in the provided dataframe: {df.index[:5]}"
        )

    # filtering the DataFrame to keep only the rows where the index can be mapped
    df = df.loc[df.index.isin(mapping_dict)]

    # renaming gene names to mapped ids using mapping dict
    df.index = df.index.map(mapping_dict)

    # TODO: check is there is another way to avoid duplicate gene names
    # sometimes different gene names have the same ensembl ID
    # for now, we just get the mean of values, but this is not ideal
    df = df.groupby(df.index).mean()

    # writing to output file
    outfile = count_file.with_name(count_file.stem + RENAMED_FILE_SUFFIX)
    df.to_csv(outfile, index=True, header=True)

    # writing mapping dict to file
    mapping_file = count_file.with_name(count_file.stem + MAPPING_FILE_SUFFIX)
    with open(mapping_file, "w") as f:
        json.dump(mapping_dict, f)


if __name__ == "__main__":
    main()
