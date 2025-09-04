#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
from parallelbar import progress_map
from Bio import Entrez
from pathlib import Path
import pandas as pd
import xmltodict
from urllib.request import urlretrieve
import tarfile
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
    before_sleep_log,
)
import yaml
from functools import partial
from multiprocessing import cpu_count
import logging

from natural_language_utils import keywords_in_fields

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ACCESSION_OUTFILE_NAME = "accessions.txt"
SPECIES_DATASETS_OUTFILE_NAME = "species_datasets.metadata.tsv"
FILTERED_DATASETS_METADATA_OUTFILE_NAME = "filtered_datasets.metadata.tsv"
REJECTED_DATASETS_METADATA_OUTFILE_NAME = "rejected_datasets.metadata.tsv"
SELECTED_DATASETS_METADATA_OUTFILE_NAME = "selected_datasets.metadata.tsv"
FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME = "selected_datasets.keywords.yaml"

ENTREZ_QUERY_MAX_RESULTS = 9999

# TODO: see how to integrate RNA-seq experiments as well
GEO_EXPERIMENT_TYPE_TO_PLATFORM = {
    "Expression profiling by array": "microarray",
    #"Expression profiling by high throughput sequencing": "rnaseq"
}

MINIML_TMPDIR = "geo_miniml"
Path(MINIML_TMPDIR).mkdir(exist_ok=True)


##################################################################
##################################################################
# EXCEPTIONS
##################################################################
##################################################################


class GeoDatasetNothingFoundError(Exception):
    pass


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
        required=True,
        help="Search GEO Datasets for this specific species"
    )
    parser.add_argument(
        "--keywords",
        type=str,
        nargs="*",
        help="Keywords to search for in datasets description",
    )
    parser.add_argument(
        "--platform",
        type=str,
        #required=True,
        help="Platform type"
    )
    parser.add_argument(
        "--exclude-accessions-in",
        dest="excluded_accessions_file",
        type=Path,
        help="Exclude accessions contained in this file",
    )
    return parser.parse_args()


@retry(
    retry=retry_if_exception_type(Exception),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def fetch_geo_datasets_for_species(species: str) -> list[dict]:
    """
    Fetch GEO datasets (GSE series) for a given species

    Args:
        species (str): Scientific name of the species (e.g. "Homo sapiens").
    """

    Entrez.email = "stableexpression@nfcore.com"
    query = f'"{species}"[Organism] AND "gse"[Entry Type]'

    # getting list of all datasets IDs for this species
    # we need possibly to perform multiple queries because the max number of returned results is capped
    nb_entries = None
    retstart = 0
    while not nb_entries or retstart < nb_entries:

        with Entrez.esearch(db="gds", term=query, retmax=ENTREZ_QUERY_MAX_RESULTS, retstart=retstart) as handle:
            record = Entrez.read(handle)

        # getting total nb of entries
        if not nb_entries:
            nb_entries = int(record["Count"])
        # setting next cursor to the next group
        retstart += ENTREZ_QUERY_MAX_RESULTS

    ids = record.get("IdList", [])
    if not ids:
        logger.warning("No GEO datasets found for your query.")
        return []

    # fetching summary info
    with Entrez.esummary(db="gds", id=",".join(ids)) as handle:
        results = Entrez.read(handle)

    # keeping only series datasets (just a double check here)
    return [ r for r in results if "GSE" in r["Accession"] ]


def exclude_unwanted_accessions(datasets: list[dict], excluded_accessions_file: Path) -> list[dict]:
    # parsing list of unwanted accessions
    with open(excluded_accessions_file) as fin:
        excluded_accessions = fin.read().splitlines()

    datasets_to_keep = []
    for dataset in datasets:
        if dataset["Accession"] not in excluded_accessions:
            datasets_to_keep.append(dataset)

    return datasets_to_keep


def format_species(species: str) -> str:
    return species.lower().replace(" ", "_")


def species_is_ok(dataset: dict, species: str) -> bool:
    accession = dataset['Accession']
    # we want datasets only specific to the species we are interested in
    parsed_species = dataset["taxon"].split("; ")
    if not parsed_species:
        logger.warning(f'Accession {accession} rejected: Could not detect species.')
        return False
    if len(parsed_species) > 1:
        logger.warning(f'Accession {accession} rejected: Found multiple species = {parsed_species}')
        return False
    if format_species(parsed_species[0]) != format_species(species):
        logger.warning(f'Accession {accession} rejected: Found wrong species = {parsed_species}')
        return False
    return True


def download_dataset_metadata(ftp_link: str, accession: str) -> Path:
    filename = f"miniml/{accession}_family.xml.tgz"
    ftp_url = ftp_link + filename
    output_file = Path(MINIML_TMPDIR) / f"{accession}.tar.gz"
    urlretrieve(ftp_url, output_file)
    return output_file


def parse_dataset_metadata(file: Path, accession: str) -> dict | None:
    with tarfile.open(file, "r:gz") as tar:
        file_to_read = f"{accession}_family.xml"
        try:
            f = tar.extractfile(file_to_read)
        except KeyError:
            file_to_read = f"{accession}_family.xml/{accession}_family.xml"
            try:
                f = tar.extractfile(file_to_read)
            except KeyError:
                return None
        if f is None:
            raise RuntimeError(f"Failed to get file: {file_to_read}")
        xml_content = f.read().decode("utf-8")
    return xmltodict.parse(xml_content)['MINiML']


def parse_characteristics(characteristics: str | dict | list, stored_characteristics: list):
    if isinstance(characteristics, str):
        stored_characteristics.append(characteristics)
    elif isinstance(characteristics, dict):
        stored_characteristics.append(characteristics["#text"])
    elif isinstance(characteristics, list):
        for c in characteristics:
            parse_characteristics(c, stored_characteristics)


def parse_interesting_metadata(dataset_metadata: dict, additional_metadata: dict) -> dict:
    sample_characteristics = []
    sample_library_strategies = []
    sample_library_sources = []
    sample_descriptions = []
    sample_titles = []
    sample_molecule_types = []

    for sample in additional_metadata["Sample"]:

        # storing description if exists
        if sample_description := sample.get("Description"):
            sample_descriptions.append(sample_description)

        # storing title if exists
        if sample_title := sample.get("Title"):
            sample_titles.append(sample_title)

            # storing molecule type if exists
            if sample_molecule_type := sample.get("Type"):
                sample_molecule_types.append(sample_molecule_type)

        # storing library strategy if exists
        if sample_library_strategy := sample.get("Library-Strategy"):
            sample_library_strategies.append(sample_library_strategy)

        # storing library source if exists
        if sample_library_source := sample.get("Library-Source"):
            sample_library_sources.append(sample_library_source)

        # parsing sample metadata
        if channels := sample.get("Channel"):
            if isinstance(channels, dict):
                channels = [channels]
            for channel in channels:
                parse_characteristics(channel["Characteristics"], sample_characteristics)

    return {
        "accession": dataset_metadata["Accession"],
        "summary": dataset_metadata["summary"],
        "title": dataset_metadata["title"],
        "overall_design": additional_metadata['Series']["Overall-Design"],
        "experiment_types": dataset_metadata["gdsType"],
        "sample_characteristics": list(set(sample_characteristics)),
        "sample_library_strategies": list(set(sample_library_strategies)),
        "sample_library_sources": list(set(sample_library_sources)),
        "sample_descriptions": list(set(sample_descriptions)),
        "sample_titles": list(set(sample_titles)),
        "sample_molecule_types": list(set(sample_molecule_types)),
    }


def format_platform_name(platform_name: str) -> str:
    return (
        platform_name
        .replace("_", "")
        .replace("-", "")
        .lower()
    )

def contains_only_rna(molecules_types: list, accession) -> bool:
    # we want only GEO series that contain only RNA molecules
    # for other series, they should be superseries contained other series that are being parsed too
    # so anyway, this would lead in duplicates
    if all([ "rna" in molecule_type.lower() for molecule_type in molecules_types]):
        return True
    logger.info(f'Accession {accession} rejected: Molecule type(s) = {molecules_types}')
    return False


def contains_proper_experiment_type(experiment_types: list, accession: str, platform: str) -> bool:
    for experiment_type in experiment_types:
        # if at least one experiment type is ok, we keep this dataset
        if GEO_EXPERIMENT_TYPE_TO_PLATFORM.get(experiment_type) == platform:
            return True
    logger.info(f'Accession {accession} rejected: Experiment type(s) = {experiment_types}')
    return False

"""
def has_proper_library_strategies(library_strategies: list, accession: str, platform: str) -> bool:
    if not library_strategies:
        logger.warning(f'No library strategies found for accession {accession}')
        # since we cannot infer, we return True
        return True

    if len(library_strategies) > 1: # multiple different platform technologies found
        logger.info(f'Multiple library strategies found for accession {accession}: {library_strategies}')
        return False

    if platform is not None:
        parsed_platform_name = library_strategies[0]
        formatted_platform_name = format_platform_name(parsed_platform_name)
        if formatted_platform_name != platform:
            logger.info(f'Accession {accession} rejected: Platform = {parsed_platform_name}')
            return False

    return True
"""


def contains_transcriptomic_source(library_sources: list, accession: str) -> bool:
    if library_sources:
        if "transcriptomic" not in library_sources:
            logger.info(f'Accession {accession} rejected: Source(s) = {library_sources}')
            return False
    return True


def parse_metadata(dataset_metadata: dict) -> dict | None:
    accession = dataset_metadata["Accession"]
    ftp_link = dataset_metadata["FTPLink"].replace("ftp://", "https://")
    downloaded_file = download_dataset_metadata(ftp_link, accession)
    additional_metadata = parse_dataset_metadata(downloaded_file, accession)

    # if we could not get additional metadata, we lack too much information to conclude
    if additional_metadata is None:
        logger.warning(f"Skipping {accession} as additional metadata is missing")
        return None

    # parsing interesting information in all available metadata
    return parse_interesting_metadata(dataset_metadata, additional_metadata)


def dataset_is_valid(metadata: dict, platform: str) -> bool:
    accession = metadata["accession"]
    # checking platform
    if not contains_proper_experiment_type(metadata["experiment_types"], accession, platform):
        return False

    # checking that library sources fit
    if not contains_transcriptomic_source(metadata["sample_library_sources"], accession):
        return False

    # checking that all molecule types are RNA
    molecules_types = metadata["sample_molecule_types"]
    if not contains_only_rna(molecules_types, accession):
        return False

    return True


def filter_metadata_with_keywords(metadata: dict, keywords: list[str]) -> dict | None:
    all_searchable_fields = (
        [metadata["summary"], metadata["title"]]
        + metadata["sample_characteristics"]
        + metadata["sample_descriptions"]
        + metadata["sample_titles"]
    )
    found_keywords = keywords_in_fields(all_searchable_fields, keywords)
    # only returning experiments if found keywords
    if found_keywords:
        metadata["found_keywords"] = list(set(found_keywords))
        logger.info(f"Found keywords: {found_keywords} in accession {metadata['accession']}")
        return metadata
    else:
        return None


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################

def main():
    args = parse_args()

    selected_accessions = []

    ncpus = cpu_count() - 1

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING GEO DATASETS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Getting datasets corresponding to species {args.species}")
    dataset_metadata_list = fetch_geo_datasets_for_species(args.species)
    logger.info(f"Found {len(dataset_metadata_list)} datasets for species {args.species}")

    #dataset_metadata_list = [d for d in dataset_metadata_list if d['Accession'] == 'GSE8203']
    #print(dataset_metadata_list[0])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXCLUDING UNWANTED ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.excluded_accessions_file:
        logger.info(f"Excluding unwanted datasets")
        dataset_metadata_list = exclude_unwanted_accessions(dataset_metadata_list, args.excluded_accessions_file)
        logger.info(f"{len(dataset_metadata_list)} datasets remaining after excluding unwanted accessions")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXCLUDING DATASETS WITH MORE THAN ONE SPECIES OR WRONG SPECIES
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logger.info(f"Excluding wrong species")
    # TODO: see how to parse data from our species from combined GEO datasets
    dataset_metadata_list = [dataset for dataset in dataset_metadata_list if species_is_ok(dataset, args.species)]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Parsing metadata")
    metadata_list = progress_map(parse_metadata, dataset_metadata_list, n_cpu=ncpus)
    metadata_list = [metadata for metadata in metadata_list if metadata is not None]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # CHECKING MOLECULE TYPE / PLATFORM TECHNOLOGIES
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info("Validating datasets")
    filtered_metadata_list = [metadata for metadata in metadata_list if dataset_is_valid(metadata, args.platform)]
    rejected_metadata_list = [metadata for metadata in metadata_list if metadata not in filtered_metadata_list]
    logger.info(f"{len(filtered_metadata_list)} datasets remaining after checking technology platform and molecule type")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FILTERING WITH KEYWORDS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.keywords:
        logger.info(f"Filtering experiments with keywords {args.keywords}")
        func = partial(filter_metadata_with_keywords, keywords=args.keywords)
        selected_metadata_list = progress_map(func, filtered_metadata_list, n_cpu=ncpus)
        selected_metadata_list = [metadata for metadata in selected_metadata_list if metadata is not None]
    else:
        selected_metadata_list = filtered_metadata_list

    if selected_metadata_list:
        logger.info(f"Kept {len(selected_metadata_list)} datasets")
        # getting accessions of selected experiments
        selected_accessions = [metadata["accession"] for metadata in selected_metadata_list]

    else:
        msg = f"Could not find experiments for species {args.species}"
        if args.keywords:
            msg += f" and keywords {args.keywords}"
        logger.warning(msg)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXPORTING DATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # exporting list of accessions
    logger.info(f"Writing accessions to {ACCESSION_OUTFILE_NAME}")
    with open(ACCESSION_OUTFILE_NAME, "w") as fout:
        fout.writelines([f"{acc}\n" for acc in selected_accessions])

    # exporting metadata
    logger.info(
        f"Writing metadata of all experiments for species {args.species} to {SPECIES_DATASETS_OUTFILE_NAME}"
    )
    df = pd.DataFrame.from_dict(dataset_metadata_list)
    df.to_csv(SPECIES_DATASETS_OUTFILE_NAME, sep="\t", index=False, header=True)

    if filtered_metadata_list:
        logger.info(f"Writing metadata of filtered datasets to {FILTERED_DATASETS_METADATA_OUTFILE_NAME}")
        df = pd.DataFrame.from_dict(filtered_metadata_list)
        df.to_csv(
            FILTERED_DATASETS_METADATA_OUTFILE_NAME,
            sep="\t",
            index=False,
            header=True,
        )

    if rejected_metadata_list:
        logger.info(f"Writing metadata of rejected datasets to {REJECTED_DATASETS_METADATA_OUTFILE_NAME}")
        df = pd.DataFrame.from_dict(rejected_metadata_list)
        df.to_csv(
            REJECTED_DATASETS_METADATA_OUTFILE_NAME,
            sep="\t",
            index=False,
            header=True,
        )

    if selected_metadata_list:
        logger.info(f"Writing metadata of selected datasets to {SELECTED_DATASETS_METADATA_OUTFILE_NAME}")
        df = pd.DataFrame.from_dict(selected_metadata_list)
        df.to_csv(
            SELECTED_DATASETS_METADATA_OUTFILE_NAME,
            sep="\t",
            index=False,
            header=True,
        )

        # exporting in YAML format too
        logger.info(
            f"Writing filtered experiments with keywords to {FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME}"
        )
        with open(FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME, "w") as fout:
            yaml.dump(selected_metadata_list, fout)


if __name__ == "__main__":
    main()

