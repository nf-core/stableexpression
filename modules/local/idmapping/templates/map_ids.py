#!/usr/bin/env python3

from gprofiler import GProfiler
import pandas as pd
from pathlib import Path

class NoIDFoundException(Exception):
    pass


gp = GProfiler(return_dataframe=True)


##################################################################
# CONSTANTS
##################################################################

RENAMED_FILE_SUFFIX = '_renamed.csv'
CHUNKSIZE = 2000

TARGET_DATABASE = 'ENSG' # Ensembl database


##################################################################
# FUNCTIONS
##################################################################

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
    splitted_species = species.lower().replace('_', ' ').split(' ')
    return splitted_species[0][0] + splitted_species[1]


def chunk_list(lst: list, chunksize: int):
    """Splits a list into chunks of a given size.

    Args:
        lst (list): The list to split.
        chunksize (int): The size of each chunk.

    Returns:
        list: A list of chunks, where each chunk is a list of len(chunksize).
    """
    return [lst[i: i + chunksize] for i in range(0, len(lst), chunksize)]


def convert_ids(gene_ids: list, species: str):

    """
    Convert a list of gene IDs to another namespace.

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
    df = gp.convert(
        organism=species,
        query=gene_ids,
        target_namespace=TARGET_DATABASE
        )

    if df.empty:
        return {}

    # keeping only rows where 'converted' is not null and only the columns of interest
    df = df.loc[df['converted'].notna(), ['incoming', 'converted']]
    df.set_index('incoming', inplace=True)
    return df.to_dict()['converted']

##################################################################
# MAIN
##################################################################




def main():

    # getting arguments
    species_name = format_species_name('$species')
    count_file = Path('$count_file')

    df = pd.read_csv(count_file, header=0, index_col=0)

    gene_ids = df.index.tolist()
    mapping_dict = {}



    chunks = chunk_list(gene_ids, chunksize=CHUNKSIZE)
    for chunk_gene_names in chunks:
        # converting to uniprot IDs / NCBI Gene IDs for all IDs comprised in this chunk
        gene_mapping = convert_ids(chunk_gene_names, species_name)
        mapping_dict.update(gene_mapping)

    if not mapping_dict: # if mapping dict is empty
        raise NoIDFoundException(
            f'No mapping found for gene names in count file {count_file.name} '
            f'and for species {species_name}! '
            f'Example of gene names: {count_file.index[:5]}')

    # filtering the DataFrame to keep only the rows where the index can be mapped
    df = df.loc[df.index.isin(mapping_dict)]

    # renaming gene names to mapped ids using mapping dict
    df.index = df.index.map(mapping_dict)

    # writing to output file
    outfile = count_file.with_name(count_file.stem + RENAMED_FILE_SUFFIX)
    df.to_csv(outfile, index=True, header=True)


if __name__ == "__main__":
    main()
