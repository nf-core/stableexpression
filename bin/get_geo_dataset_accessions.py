#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import random
import tarfile
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from urllib.request import urlretrieve

import pandas as pd
import requests
import xmltodict
from Bio import Entrez
from natural_language_utils import keywords_in_fields
from requests.exceptions import ConnectionError, HTTPError
from resource_management import set_max_memory
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# set a custom writable directory before any Entrez operations
# mandatory for running the script in an apptainer container
# Entrez.Parser.Parser.directory("/tmp/biopython")

ALLOWED_PLATFORMS = ["rnaseq", "microarray"]

ACCESSION_OUTFILE_NAME = "accessions.txt"
SPECIES_DATASETS_OUTFILE_NAME = "geo_all_datasets.metadata.tsv"
REJECTED_DATASETS_OUTFILE_NAME = "geo_rejected_datasets.metadata.tsv"
# WRONG_SPECS_DATASETS_METADATA_OUTFILE_NAME = "geo_wrong_platform_moltype_datasets.metadata.tsv"
# WRONG_KEYWORDS_DATASETS_METADATA_OUTFILE_NAME = "geo_wrong_keywords_datasets.metadata.tsv"
# PLATFORM_NOT_AVAILABLE_DATASETS_METADATA_OUTFILE_NAME = "platform_not_available_datasets.metadata.tsv"
# GENE_ID_MAPPING_ISSUES_DATASETS_METADATA_OUTFILE_NAME = "gene_id_mapping_issues_datasets.metadata.tsv"
SELECTED_DATASETS_OUTFILE_NAME = "geo_selected_datasets.metadata.tsv"

ENTREZ_QUERY_MAX_RESULTS = 9999
ENTREZ_EMAIL = "stableexpression@nfcore.com"
PLATFORM_METADATA_CHUNKSIZE = 2000

NCBI_API_BASE_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc={accession}"
)
STOP_RETRY_AFTER_DELAY = 600

NB_PROBE_IDS_TO_PARSE = 1000
NB_PROBE_IDS_TO_SAMPLE = 10

SUPERSERIES_SUMMARY = "This SuperSeries is composed of the SubSeries listed below."

ALLOWED_LIBRARY_SOURCES = ["transcriptomic", "RNA"]
ALLOWED_MOLECULE_TYPES = ["RNA", "SRA"]

GEO_EXPERIMENT_TYPE_TO_PLATFORM = {
    "Expression profiling by array": "microarray",
    "Expression profiling by high throughput sequencing": "rnaseq",
}

MINIML_TMPDIR = "geo_miniml"
PLATFORM_SOFT_TMPDIR = "geo_platform_soft"
Path(MINIML_TMPDIR).mkdir(exist_ok=True)
Path(PLATFORM_SOFT_TMPDIR).mkdir(exist_ok=True)


##################################################################
##################################################################
# EXCEPTIONS
##################################################################
##################################################################


class GeoDatasetNothingFoundError(Exception):
    pass


class GeoPlatformDataTableNotFound(Exception):
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
        help="Search GEO Datasets for this specific species",
    )
    parser.add_argument(
        "--keywords",
        type=str,
        nargs="*",
        help="Keywords to search for in datasets description",
    )
    parser.add_argument(
        "--platform", type=str, help="Platform type", choices=ALLOWED_PLATFORMS
    )
    parser.add_argument(
        "--exclude-accessions-in",
        dest="excluded_accessions_file",
        type=Path,
        help="Exclude accessions contained in this file",
    )
    parser.add_argument(
        "--random-sampling-size",
        dest="random_sampling_size",
        type=int,
        help="Random sampling size",
    )
    parser.add_argument(
        "--random-sampling-seed",
        dest="random_sampling_seed",
        type=int,
        help="Random sampling seed",
    )
    parser.add_argument(
        "--cpus", type=str, dest="nb_cpus", required=True, help="Number of CPUs"
    )
    parser.add_argument(
        "--memory", type=str, dest="memory", required=True, help="Memory in GB"
    )
    parser.add_argument(
        "--accessions",
        type=str,
        help="[For dev purposes / testing: provide directly accessions (separated by commas) and try to get their metadata]",
    )
    return parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUERIES TO ENTREZ
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
    retry_error_callback=(lambda _: {}),
)
def send_request_to_entrez_esearch(query: str) -> dict:
    Entrez.email = ENTREZ_EMAIL
    with Entrez.esearch(
        db="gds", term=query, retmax=ENTREZ_QUERY_MAX_RESULTS
    ) as handle:
        return Entrez.read(handle)


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
    retry_error_callback=(lambda _: []),
)
def send_request_to_entrez_esummary(ids: list[str]) -> list[dict]:
    Entrez.email = ENTREZ_EMAIL
    ids_str = ",".join(ids)
    with Entrez.esummary(
        db="gds", id=ids_str, retmax=ENTREZ_QUERY_MAX_RESULTS
    ) as handle:
        return Entrez.read(handle)


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
    retry_error_callback=(lambda _: None),
)
def send_request_to_ncbi_api(accession: str) -> requests.Response | None:
    url = NCBI_API_BASE_URL.format(accession=accession)
    server_error = False
    response = None

    try:
        response = requests.get(url, stream=True)
    except requests.exceptions.ConnectionError:
        server_error = True
    else:
        try:
            response.raise_for_status()
        except (HTTPError, ConnectionError) as err:
            if str(response.status_code).startswith("5"):  # error 500 -> 509
                server_error = True
                raise err
            else:
                logger.error(
                    f"Error {response.status_code} while sending request to NCBI: {err}"
                )
                raise err

    # if we get connection issues or 500 -> 509 server errors
    # we stop immediately for this accession (return None)
    if server_error:
        logger.critical(
            f"Server error while sending request to NCBI for accession {accession}"
        )

    return response


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
    retry_error_callback=(lambda _: None),
)
def download_file_at_url(url: str, output_file: Path):
    urlretrieve(url, output_file)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GEO DATASETS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def chunk_list(lst: list, chunksize: int) -> list:
    """Splits a list into chunks of a given size.

    Args:
        lst (list): The list to split.
        chunksize (int): The size of each chunk.

    Returns:
        list: A list of chunks, where each chunk is a list of len(chunksize).
    """
    return [lst[i : i + chunksize] for i in range(0, len(lst), chunksize)]


def fetch_geo_datasets_for_species(species: str) -> list[dict]:
    """
    Fetch GEO datasets (GSE series) for a given species

    Args:
        species (str): Scientific name of the species (e.g. "Homo sapiens").
    """
    dataset_types = [
        f'"{experiment_type}"[DataSet Type]'
        for experiment_type in GEO_EXPERIMENT_TYPE_TO_PLATFORM
    ]
    formatted_dataset_type = "(" + " OR ".join(dataset_types) + ")"

    query = f'"{species}"[Organism] AND "gse"[Entry Type] AND {formatted_dataset_type}'
    logger.info(f"Fetching GEO datasets with query: {query}")

    # getting list of all datasets IDs for this species
    # we need possibly to perform multiple queries because the max number of returned results is capped
    nb_entries = None
    retstart = 0
    while not nb_entries or retstart < nb_entries:
        record = send_request_to_entrez_esearch(query)

        if not record:
            logger.warning(f"Failed to query Entrey Esearch with query: {query}")
            return []

        # getting total nb of entries
        if not nb_entries:
            nb_entries = int(record["Count"])

            # if there is no entry for this species
            if nb_entries == 0:
                logger.info(f"No entries found for query: {query}")
                return []

        # setting next cursor to the next group
        retstart += ENTREZ_QUERY_MAX_RESULTS

    ids = record.get("IdList", [])
    if not ids:
        logger.warning("No GEO datasets found for your query.")
        return []

    # fetching summary info
    results = send_request_to_entrez_esummary(ids)

    # keeping only series datasets (just a double check here)
    # and removing superseries (they are just containers of series that are also contained here)
    return [
        r
        for r in results
        if "GSE" in r["Accession"] and r["summary"] != SUPERSERIES_SUMMARY
    ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FORMATTING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def format_species(species: str) -> str:
    return "_".join(species.lower().split(" "))


def format_platform_name(platform_name: str) -> str:
    return platform_name.replace("_", "").replace("-", "").lower()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET METADATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def download_dataset_metadata(ftp_link: str, accession: str) -> Path | None:
    filename = f"miniml/{accession}_family.xml.tgz"
    ftp_url = ftp_link + filename
    output_file = Path(MINIML_TMPDIR) / f"{accession}.tar.gz"
    download_file_at_url(ftp_url, output_file)
    if output_file.exists():
        return output_file
    else:
        logger.error(f"Failed to download dataset metadata for accession: {accession}")
        return None


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
            logger.warning(f"Failed to get file: {file_to_read}")
            return None

        try:
            xml_content = f.read().decode("utf-8")
        except UnicodeDecodeError:
            logger.warning(f"Failed to decode file: {file_to_read}")
            return None

    return xmltodict.parse(xml_content)["MINiML"]


def parse_characteristics(
    characteristics: str | dict | list, stored_characteristics: list
):
    if isinstance(characteristics, str):
        stored_characteristics.append(characteristics)
    elif isinstance(characteristics, dict):
        if "#text" in characteristics:
            stored_characteristics.append(characteristics["#text"])
    elif isinstance(characteristics, list):
        for c in characteristics:
            parse_characteristics(c, stored_characteristics)


def parse_interesting_metadata(
    dataset_metadata: dict, additional_metadata: dict
) -> dict:
    """
    Parses interesting metadata from a dataset metadata dictionary and additional metadata dictionary.

    Args:
        dataset_metadata (dict): The dataset metadata dictionary.
        additional_metadata (dict): The additional metadata dictionary.

    Returns:
        dict: The parsed interesting metadata dictionary.
    """
    sample_characteristics = []
    sample_library_strategies = []
    sample_library_sources = []
    sample_descriptions = []
    sample_titles = []
    sample_molecule_types = []

    platform_accessions = [
        "GPL" + gpl_id for gpl_id in dataset_metadata["GPL"].split(";")
    ]

    experiment_types = dataset_metadata["gdsType"]
    experiment_types = (
        experiment_types if isinstance(experiment_types, list) else [experiment_types]
    )

    # if additional metadata have sample information
    if "Sample" in additional_metadata:
        # change to list if it's a single dictionary
        if isinstance(additional_metadata["Sample"], dict):
            additional_metadata["Sample"] = [additional_metadata["Sample"]]

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
                    parse_characteristics(
                        channel["Characteristics"], sample_characteristics
                    )

    return {
        "accession": dataset_metadata["Accession"],
        "taxon": dataset_metadata["taxon"],
        "platform_accessions": platform_accessions,
        "summary": dataset_metadata["summary"],
        "title": dataset_metadata["title"],
        "overall_design": additional_metadata["Series"]["Overall-Design"],
        "experiment_types": experiment_types,
        "sample_characteristics": list(set(sample_characteristics)),
        "sample_library_strategies": list(set(sample_library_strategies)),
        "sample_library_sources": list(set(sample_library_sources)),
        "sample_descriptions": list(set(sample_descriptions)),
        "sample_titles": list(set(sample_titles)),
        "sample_molecule_types": list(set(sample_molecule_types)),
    }


def fetch_dataset_metadata(dataset_metadata: dict) -> dict | None:
    """
    Parses metadata from a dataset metadata dictionary.

    Args:
        dataset_metadata (dict): The dataset metadata dictionary.

    Returns:
        dict | None: The parsed metadata dictionary or None if the metadata is missing.
    """
    accession = dataset_metadata["Accession"]
    ftp_link = dataset_metadata["FTPLink"].replace("ftp://", "https://")
    downloaded_file = download_dataset_metadata(ftp_link, accession)
    if downloaded_file is None:
        logger.warning(f"Skipping {accession} as metadata download failed")
        return None

    additional_metadata = parse_dataset_metadata(downloaded_file, accession)

    # if we could not get additional metadata, we lack too much information to conclude
    if additional_metadata is None:
        logger.warning(f"Skipping {accession} as additional metadata is missing")
        return None

    # parsing interesting information in all available metadata
    return parse_interesting_metadata(dataset_metadata, additional_metadata)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# METADATA TESTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def exclude_unwanted_accessions(
    datasets: list[dict], excluded_accessions: list[str]
) -> tuple[list[dict], list[dict]]:
    datasets_to_keep = []
    excluded_datasets = []
    for dataset in datasets:
        if dataset["accession"] in excluded_accessions:
            excluded_datasets.append(dataset)
        else:
            datasets_to_keep.append(dataset)
    return datasets_to_keep, excluded_datasets


def check_species_issues(parsed_species_list: list, species: str) -> str | None:
    # trying to find our species in the list of species parsed
    for parsed_species in parsed_species_list:
        if format_species(parsed_species) == format_species(species):
            return None
    return f"PARSED SPECIES: {parsed_species_list}"


def check_molecule_type_issues(molecules_types: list) -> str | None:
    # we want only GEO series that contain only RNA molecules
    # for other series, they should be superseries contained other series that are being parsed too
    # so anyway, this would lead in duplicates
    if any(
        [
            molecule_type.upper() in ALLOWED_MOLECULE_TYPES
            for molecule_type in molecules_types
        ]
    ):
        return None
    return f"MOLECULE TYPES: {molecules_types}"


def check_experiment_type_issues(experiment_types: list, platform: str) -> str | None:
    for experiment_type in experiment_types:
        # if at least one experiment type is ok, we keep this dataset
        if GEO_EXPERIMENT_TYPE_TO_PLATFORM.get(experiment_type) == platform:
            return None
    return f"EXPERIMENT TYPES: {experiment_types}"


def check_source_issues(library_sources: list) -> str | None:
    # if we have no data about library sources, we just cannot infer
    if not library_sources:
        return None
    if any(
        library_source in ALLOWED_LIBRARY_SOURCES for library_source in library_sources
    ):
        return None
    return f"LIBRARY SOURCES: {library_sources}"


def search_keywords(dataset: dict, keywords: list[str]) -> tuple[list, str | None]:
    accession = dataset["accession"]
    all_searchable_fields = (
        [dataset["summary"], dataset["title"]]
        + dataset["sample_characteristics"]
        + dataset["sample_descriptions"]
        + dataset["sample_titles"]
    )
    found_keywords = keywords_in_fields(all_searchable_fields, keywords)
    # only returning experiments if found keywords
    if found_keywords:
        dataset["found_keywords"] = list(set(found_keywords))
        logger.info(f"Found keywords: {found_keywords} in accession {accession}")
        return found_keywords, None
    else:
        return [], "NO KEYWORDS_FOUND"


def check_dataset(
    dataset: dict, species: str, platform: str | None, keywords: list[str] | None
) -> tuple[list, dict]:
    accession = dataset["accession"]
    parsed_species_list = dataset["taxon"].split("; ")
    experiment_types = dataset["experiment_types"]
    library_sources = dataset["sample_library_sources"]
    molecules_types = dataset["sample_molecule_types"]

    issues = []

    # checking species
    if issue := check_species_issues(parsed_species_list, species):
        issues.append(issue)

    # checking platform
    if platform is not None:
        if issue := check_experiment_type_issues(experiment_types, platform):
            issues.append(issue)

    # checking that library sources fit
    if issue := check_source_issues(library_sources):
        issues.append(issue)

    # checking that all molecule types are RNA
    if issue := check_molecule_type_issues(molecules_types):
        issues.append(issue)

    found_keywords = []
    if keywords:
        found_keywords, keyword_issue = search_keywords(dataset, keywords)
        if keyword_issue:
            issues.append(keyword_issue)

    if issues:
        rejection_dict = {accession: issues}
    else:
        rejection_dict = {}

    return found_keywords, rejection_dict


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GEO PLATFORMS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def fetch_geo_platform_metadata(datasets: list[dict]) -> dict:
    """
    Fetch data for a GEO platform

    Args:
        platform_accession (str): accession of the platform
    """
    # unique list of platform accessions
    platform_accessions = list(
        set(
            [
                platform_accession
                for dataset in datasets
                for platform_accession in dataset["platform_accessions"]
            ]
        )
    )
    # formating query
    formatted_platform_accessions = [
        f'"{platform_accession}"[GEO Accession]'
        for platform_accession in platform_accessions
    ]
    platform_accessions_str = " OR ".join(formatted_platform_accessions)
    query = f'({platform_accessions_str}) AND "gpl"[Entry Type] '

    record = send_request_to_entrez_esearch(query=query)

    ids = record.get("IdList", [])
    if not ids:
        logger.warning(f"No GEO platform found for accessions {platform_accessions}.")
        return {}

    # fetching summary info
    # one single request to NCBI for all platform accessions
    platform_metadatas = send_request_to_entrez_esummary(ids)
    # return dict associating dataset accessions with platform metadata
    return {
        platform_metadata["Accession"]: platform_metadata
        for platform_metadata in platform_metadatas
    }


def check_dataset_platforms(
    dataset: dict, accession_to_platform_metadata: dict, species: str
) -> dict:
    accession = dataset["accession"]
    platform_accessions = dataset["platform_accessions"]

    if not platform_accessions:
        return {accession: "NO PLATFORM ACCESSIONS"}

    platforms_metadata = [
        accession_to_platform_metadata[platform_accession]
        for platform_accession in dataset["platform_accessions"]
    ]

    # getting list of platform taxon
    platforms_taxons = []
    for metadata in platforms_metadata:
        if metadata.get("taxon") is not None:
            platforms_taxons += metadata.get("taxon").split("; ")
    platforms_taxons = list(set(platforms_taxons))

    if not platforms_taxons:
        return {accession: "NO PLATFORM TAXON"}

    # checking if at least one of the platform accession is the good one
    # sample will be further filtered during download (download_geo_data.R)
    if not any(
        format_species(species) == format_species(taxon) for taxon in platforms_taxons
    ):
        return {accession: f"TAXON MISMATCH: {platforms_taxons}"}

    return {}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RANDOM SAMPLING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def sample_experiments_randomly(
    experiments: list[dict], sampling_size: int, seed: int
) -> list[str]:
    random.seed(seed)
    sampled_experiments = []

    total_nb_samples = 0
    experiments_left = list(experiments)
    while experiments_left and total_nb_samples <= sampling_size:
        # if the min number of samples is greater than the remaining space left, we get out of the loop
        experiments_left_nb_samples = [exp["nb_samples"] for exp in experiments_left]
        min_nb_samples = min(experiments_left_nb_samples)
        if min_nb_samples > sampling_size - total_nb_samples:
            break

        found_experiment = False
        test_total_nb_samples = int(total_nb_samples)
        not_chosen_yet = list(experiments_left)
        while not_chosen_yet and not found_experiment:
            experiment = random.choice(not_chosen_yet)
            not_chosen_yet.remove(experiment)
            test_total_nb_samples = total_nb_samples + experiment["nb_samples"]
            if test_total_nb_samples <= sampling_size:
                found_experiment = True

        # if the last one was not good, it means we reached the limit of samples we can take
        if not found_experiment:
            break
        else:
            total_nb_samples = test_total_nb_samples
            experiments_left.remove(experiment)
            sampled_experiments.append(experiment)

    return [exp["accession"] for exp in sampled_experiments]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def sort_if_list(x):
    if isinstance(x, list):
        return sorted(x)
    else:
        return x


def export_dataset_metadatas(
    datasets: list[dict], output_file: str, clean_columns: bool = True
):
    if datasets:
        df = pd.DataFrame.from_dict(datasets)
        # all dataframe contain the column "accession"
        # sorting by accessions to ensure that outputs are reproducible
        df.sort_values(by="accession", inplace=True)
        for col in df.columns:
            df[col] = df[col].apply(sort_if_list)
        # cleaning columns so that MultiQC can parse them
        if clean_columns:
            for col in df.columns:
                df[col] = df[col].astype(str).str.replace("\n", "")
                df[col] = df[col].astype(str).str.replace("\t", "")
        df.to_csv(
            output_file,
            sep="\t",
            index=False,
            header=True,
        )


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    args = parse_args()

    set_max_memory(args.memory)

    random_sampling_size = args.random_sampling_size

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING GEO DATASETS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Getting datasets corresponding to species {args.species}")
    datasets = fetch_geo_datasets_for_species(args.species)
    logger.info(f"Found {len(datasets)} datasets for species {args.species}")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FOR DEV PURPOSES / TESTING: RESTRICT TO SPECIFIC ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.accessions:
        logger.info(f"Keeping only accessions {args.accessions}")
        dev_accessions = args.accessions.split(",")
        datasets = [d for d in datasets if d["Accession"] in dev_accessions]
        logger.info(f"Kept {len(datasets)} datasets for dev / testing purposes")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING DATASET METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Parsing metadata for {len(datasets)} datasets")
    augmented_datasets = []
    with (
        Pool(processes=args.nb_cpus) as p,
        tqdm(total=len(datasets)) as pbar,
    ):
        for result in p.imap_unordered(fetch_dataset_metadata, datasets):
            pbar.update()
            pbar.refresh()
            if result is None:
                continue
            augmented_datasets.append(result)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # VALIDATING DATASETS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Validating {len(augmented_datasets)} datasets")
    checked_datasets = []
    rejection_dict = {}
    for dataset in tqdm(augmented_datasets):
        found_keywords, issue_dict = check_dataset(
            dataset, args.species, args.platform, args.keywords
        )
        if issue_dict:
            rejection_dict |= issue_dict
        else:
            if found_keywords:
                dataset["found_keywords"] = found_keywords
            checked_datasets.append(dataset)

    logger.info(f"Validated {len(checked_datasets)} datasets")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXCLUDING UNWANTED ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # we exclude unwanted accessions only now
    # because we want to get the metadata of the excluded datasets
    # in order to adjust the random sampling size
    if args.excluded_accessions_file:
        # parsing list of accessions which were already fetched from Expression Atlas
        with open(args.excluded_accessions_file) as fin:
            excluded_accessions = fin.read().splitlines()
        logger.info("Excluding unwanted datasets")
        checked_datasets, excluded_datasets = exclude_unwanted_accessions(
            checked_datasets, excluded_accessions
        )
        logger.info(
            f"{len(checked_datasets)} datasets remaining after excluding unwanted accessions"
        )

        # adjusting random sampling size by substracting the number of excluded accessions
        if random_sampling_size:
            total_nb_excluded_samples = sum(
                [len(dataset["sample_titles"]) for dataset in excluded_datasets]
            )
            logger.info(
                f"Subtracting {total_nb_excluded_samples} samples from random sampling size"
            )
            random_sampling_size -= total_nb_excluded_samples
            # keeping it positive (just in case)
            if random_sampling_size < 0:
                logger.warning(
                    f"Random sampling size is negative ({random_sampling_size}), setting it to 0"
                )
                random_sampling_size = 0

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GETTING METADATA OF SEQUENCING PLATFORMS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info("Getting platform metadata")
    # making chunks to group requests to NCBI GEO
    checked_datasets_chunks = chunk_list(checked_datasets, PLATFORM_METADATA_CHUNKSIZE)
    # resetting selecting datasets
    accession_to_platform_metadata = {}
    for selected_datasets_chunk in tqdm(checked_datasets_chunks):
        accession_to_platform_metadata |= fetch_geo_platform_metadata(
            selected_datasets_chunk
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # VALIDATING EACH PLATFORM SEPARATELY, DATASET BY DATASET
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Checking each platform for {len(checked_datasets)} datasets")
    func = partial(
        check_dataset_platforms,
        accession_to_platform_metadata=accession_to_platform_metadata,
        species=args.species,
    )
    selected_datasets = []
    # resetting selecting datasets
    for dataset in tqdm(checked_datasets):
        accession = dataset["accession"]
        issue_dict = func(dataset)
        if issue_dict:
            if accession in rejection_dict:  # should not happen but in case
                rejection_dict[accession] += issue_dict[accession]
            else:
                rejection_dict |= issue_dict
        else:
            selected_datasets.append(dataset)

    if rejection_dict:
        logger.warning(f"{len(rejection_dict)} datasets rejected")
        logger.warning(f"Reasons for rejection: {rejection_dict}")

    selected_accessions = sorted(
        [dataset["accession"] for dataset in selected_datasets]
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # RANDOM SAMPLING
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if random_sampling_size is not None and args.random_sampling_seed is not None:
        selected_accession_to_nb_samples = [
            {
                "accession": dataset["accession"],
                "nb_samples": len(dataset["sample_titles"]),
            }
            for dataset in selected_datasets
        ]

        nb_samples_df = pd.DataFrame.from_dict(selected_accession_to_nb_samples)
        nb_samples_df.to_csv("selected_accession_to_nb_samples.csv", index=False)

        logger.info("Sampling experiments randomly")
        selected_accessions = sample_experiments_randomly(
            selected_accession_to_nb_samples,
            random_sampling_size,
            args.random_sampling_seed,
        )
        logger.info(
            f"Kept {len(selected_accessions)} experiments after random sampling"
        )
        selected_datasets = [
            dataset
            for dataset in selected_datasets
            if dataset["accession"] in selected_accessions
        ]
    else:
        logger.info(
            f"No random sampling requested. Kept {len(selected_datasets)} datasets"
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXPORTING ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # sorting accessions to ensure that outputs are reproducible
    selected_accessions = sorted(selected_accessions)
    with open(ACCESSION_OUTFILE_NAME, "w") as fout:
        fout.write("\n".join(selected_accessions))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXPORTING DATASETS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    export_dataset_metadatas(augmented_datasets, SPECIES_DATASETS_OUTFILE_NAME)
    export_dataset_metadatas(selected_datasets, SELECTED_DATASETS_OUTFILE_NAME)

    rejected_datasets = [
        {"accession": accession, "reason": reason}
        for accession, reason in rejection_dict.items()
    ]
    export_dataset_metadatas(
        rejected_datasets, REJECTED_DATASETS_OUTFILE_NAME, clean_columns=False
    )


if __name__ == "__main__":
    main()
