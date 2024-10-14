#!/usr/bin/env python3
from Bio import Entrez
import pandas as pd
from functools import partial
from pathlib import Path
from retry import retry
from urllib.error import HTTPError


class EsearchHTTPError(Exception):
    pass

##################################################################
##################################################################
# CONSTANTS
##################################################################
##################################################################

EMAIL_DOMAIN = "@nfcore-sampleexpression.com"

# not more than 10000! Entrez.esearch's retmax parameters cannot go over 10000
# https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
CHUNKSIZE = 500

RENAMED_FILE_SUFFIX = '_renamed.csv'


##################################################################
##################################################################
# FUNCTIONS
##################################################################
##################################################################

def chunk_list(lst: list, chunksize: int):
    """Splits a list into chunks of a given size.

    Args:
        lst (list): The list to split.
        chunksize (int): The size of each chunk.

    Returns:
        list: A list of chunks, where each chunk is a list of len(chunksize).
    """
    return [lst[i: i + chunksize] for i in range(0, len(lst), chunksize)]


def format_species_name(species: str):
    return species.replace('_', ' ').capitalize()



@retry(EsearchHTTPError, tries=5, delay=2, backoff=2)
def get_esearch_result(search_term: str, retmax: int):
    """
    Perform an Entrez esearch query and return the search results.

    The function is decorated with a retry decorator. If the NCBI server returns
    a 429 "Too Many Requests" error, the function will wait and retry up to 5
    times with an exponential backoff delay.

    Args:
        search_term (str): The search term to query Entrez with.
        retmax (int): The maximum number of records to return.

    Returns:
        dict: The search results as a dictionary.
    """
    try:
        with Entrez.esearch(db="gene", term=search_term, retmode="xml", retmax=retmax) as handle:
            search_results = Entrez.read(handle)
    except HTTPError as err:
        print(dir(HTTPError))
        if err.status_code == 429:
            raise EsearchHTTPError
    return search_results


def get_ncbi_gene_id_mapping(gene_names: list[str], species_name: str):
    """Maps a list of gene names to their corresponding NCBI Gene IDs.

    This function first constructs an Entrez search query by combining the gene
    names with "OR" and adding the species name as a filter. It then performs the
    search and fetches the detailed information for the gene IDs. The function
    returns a dictionary mapping the gene names to their corresponding NCBI Gene
    IDs.

    Args:
        species_name (str): The name of the species to filter the results with.
        gene_names (list[str]): A list of gene names to map to their NCBI Gene IDs.

    Returns:
        dict: A dictionary mapping the gene names to their corresponding NCBI Gene
            IDs.
    """
    # making search query string
    search_term = " OR ".join([f"{gene}[Gene Name]" for gene in gene_names])  # Combine gene names with OR
    search_term += f" AND {species_name}[Organism]"

    search_results = get_esearch_result(search_term, retmax=len(gene_names))

    # Step 2: Fetch detailed information about the gene IDs
    gene_ids = search_results["IdList"]

    # Fetch detailed information for each gene ID using efetch or esummary
    with Entrez.efetch(db="gene", id=",".join(gene_ids), retmode="xml") as handle:
        gene_records = Entrez.read(handle)

    # just in case
    if len(gene_records) < len(gene_names):
        raise RuntimeError('Did not fetch as many records as gene names!')

    # Step 3: Map the gene names to their corresponding IDs by retrieving the official gene symbol
    gene_mapping = {}
    for record in gene_records:

        gene_ref_string_values = extract_strings(record['Entrezgene_gene']['Gene-ref'])

        gene_name = None
        for string in gene_ref_string_values:
            if string in gene_names:
                gene_name = string
        if gene_name is None:
            raise ValueError(f'Could not find gene name in {gene_ref_string_values}')

        gene_id = record['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
        gene_mapping[gene_name] = gene_id

    return gene_mapping


def extract_strings(d):
    """
    Extracts all strings from a nested dictionary or list.

    Args:
        d: The dictionary or list to extract strings from.

    Returns:
        A list of strings.
    """
    # List to collect all string values
    strings = []

    # Recursive function to traverse the dictionary or list
    def recursive_extractor(value):
        if isinstance(value, dict):  # If it's a dictionary, iterate over its values
            for v in value.values():
                recursive_extractor(v)
        elif isinstance(value, list):  # If it's a list, iterate over its elements
            for item in value:
                recursive_extractor(item)
        elif isinstance(value, str):  # If it's a string, collect it
            strings.append(value)

    # Start the recursion with the input dictionary
    recursive_extractor(d)

    return strings


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################

def main():

    # getting arguments
    species_name = format_species_name('$species')
    count_file = Path('$count_file')

    # (required by NCBI) using fake email here
    Entrez.email = count_file.stem + EMAIL_DOMAIN

    df = pd.read_csv(count_file, header=0, index_col=0)

    gene_ids = df.index.tolist()
    chunks = chunk_list(gene_ids, chunksize=CHUNKSIZE)

    mapping_dict = {}
    func = partial(get_ncbi_gene_id_mapping, species_name=species_name)
    for chunk_gene_names in chunks:
        # converting to uniprot IDs for all IDs comprised in this chunk
        gene_mapping = func(chunk_gene_names)
        mapping_dict.update(gene_mapping)

    # renaming gene names to ncbi gene ids using mapping dict
    df.index = df.index.map(lambda x: mapping_dict.get(x, x))

    # writing to output file
    outfile = count_file.with_suffix(RENAMED_FILE_SUFFIX)
    df.to_csv(outfile, index=True, header=True)


if __name__ == "__main__":
    main()
