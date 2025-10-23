#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
from tqdm import tqdm
from multiprocessing import Pool
from Bio import Entrez
from pathlib import Path
from random import sample
import re
import requests
import pandas as pd
import xmltodict
from urllib.request import urlretrieve
import tarfile
from tenacity import (
    retry,
    stop_after_delay,
    wait_exponential,
    before_sleep_log,
)
import yaml
from functools import partial
import logging
from requests.exceptions import HTTPError, ConnectionError

from natural_language_utils import keywords_in_fields
from gprofiler_utils import convert_ids, chunk_list

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# set a custom writable directory before any Entrez operations
# mandatory for running the script in an apptainer container
# Entrez.Parser.Parser.directory("/tmp/biopython")

ACCESSION_OUTFILE_NAME = "accessions.txt"
SPECIES_DATASETS_OUTFILE_NAME = "species_datasets.metadata.tsv"
FILTERED_DATASETS_METADATA_OUTFILE_NAME = "filtered_datasets.metadata.tsv"
REJECTED_DATASETS_METADATA_OUTFILE_NAME = "rejected_datasets.metadata.tsv"
SELECTED_DATASETS_METADATA_OUTFILE_NAME = "selected_datasets.metadata.tsv"
FINAL_DATASETS_METADATA_OUTFILE_NAME = "final_datasets.metadata.tsv"
FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME = "selected_datasets.keywords.yaml"

ENTREZ_QUERY_MAX_RESULTS = 9999
ENTREZ_EMAIL = "stableexpression@nfcore.com"
ENTREZ_CHUNKSIZE = 2000

NCBI_API_BASE_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc={accession}"
)

STOP_RETRY_AFTER_DELAY = 600

NB_PROBE_IDS_TO_PARSE = 1000
NB_PROBE_IDS_TO_SAMPLE = 10

ALLOWED_LIBRARY_SOURCES = ["transcriptomic", "RNA"]

# TODO: see how to integrate RNA-seq experiments as well
GEO_EXPERIMENT_TYPE_TO_PLATFORM = {
    "Expression profiling by array": "microarray",
    # "Expression profiling by high throughput sequencing": "rnaseq"
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
    parser.add_argument("--platform", type=str, help="Platform type")
    parser.add_argument(
        "--exclude-accessions-in",
        dest="excluded_accessions_file",
        type=Path,
        help="Exclude accessions contained in this file",
    )
    parser.add_argument(
        "--cpus",
        dest="nb_cpus",
        type=int,
        required=True,
        help="Number of CPUs to use",
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


def fetch_geo_datasets_for_species(species: str) -> list[dict]:
    """
    Fetch GEO datasets (GSE series) for a given species

    Args:
        species (str): Scientific name of the species (e.g. "Homo sapiens").
    """
    query = f'"{species}"[Organism] AND "gse"[Entry Type] AND "expression profiling by array"[DataSet Type]'
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
    return [r for r in results if "GSE" in r["Accession"]]


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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GEO PLATFORMS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def fetch_geo_platform_data(platform_accessions: list[str]) -> dict:
    """
    Fetch data for a GEO platform

    Args:
        platform_accession (str): accession of the platform
    """
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
        return []

    # fetching summary info
    results = send_request_to_entrez_esummary(ids)
    return results


def get_platform_metadata(selected_metadata_chunk_list: list[dict]) -> list[dict]:
    # unique list of platform accessions
    platform_accessions = list(
        set(
            [
                platform_accession
                for metadata in selected_metadata_chunk_list
                for platform_accession in metadata["platform_accessions"]
            ]
        )
    )
    # one single request to NCBI for all platform accessions
    # we extract the platform accessions to allow better parsing afterwards
    pltf_acc_to_pltf_metadata = {
        platform_metadata["Accession"]: platform_metadata
        for platform_metadata in fetch_geo_platform_data(platform_accessions)
    }

    # adding the platform metadata to the corresponding metadata
    augmented_metadata_list = []
    for metadata in selected_metadata_chunk_list:
        metadata["platform_metadata"] = []
        for platform_accession in metadata["platform_accessions"]:
            # augmenting metadata with platform metadata
            # filtering out cases where the platform metadata is not available
            if platform_accession in pltf_acc_to_pltf_metadata:
                metadata["platform_metadata"].append(
                    pltf_acc_to_pltf_metadata[platform_accession]
                )
        augmented_metadata_list.append(metadata)

    return augmented_metadata_list


"""
def download_platform_datatable(ftp_link: str, platform_accession: str) -> Path | None:
    filename = f"soft/{platform_accession}_family.soft.gz"
    ftp_url = ftp_link + filename
    output_file = Path(PLATFORM_SOFT_TMPDIR) / f"{platform_accession}.gz"
    download_file_at_url(ftp_url, output_file)
    return output_file
"""


def get_platform_probe_id_samples(platform_accession: str) -> list[str]:
    response = send_request_to_ncbi_api(platform_accession)
    if response is None:
        return []

    header_found = False
    probe_ids = []
    counter = 0
    for line in response.iter_lines(decode_unicode=True):
        if counter >= NB_PROBE_IDS_TO_PARSE:
            break
        if line:
            # removing HTML patterns
            line = re.sub("<[^<]+?>", "", line).strip()
            line = re.sub(r"</?strong>", "", line)
            # first things first: try to get the header
            if not header_found:
                if line.startswith("ID"):
                    header_found = True
                    continue
            else:
                # once the header was gotten, all successive lines are the data
                probe_id = line.split("\t")[0]
                probe_ids.append(probe_id)
                counter += 1

    # return a random sample of probe IDs
    nb_samples = min(len(probe_ids), NB_PROBE_IDS_TO_SAMPLE)
    return sample(probe_ids, nb_samples)


def probe_ids_can_be_converted(
    dataset_metadata: dict, species: str
) -> tuple[dict, bool]:
    platform_dict_list = dataset_metadata["platform_metadata"]
    all_probe_ids = []

    for platform_dict in platform_dict_list:
        # looping until we find data for our species
        if format_species(platform_dict["taxon"]) != format_species(species):
            continue
        # getting a sample of the first probe ids
        sampled_probe_ids = get_platform_probe_id_samples(platform_dict["Accession"])

        # if we could not get any probe ids for a platform, we won't use this dataset
        if not sampled_probe_ids:
            return dataset_metadata, False

        all_probe_ids += sampled_probe_ids

    # try to convert ids
    mapping_dict, _ = convert_ids(all_probe_ids, species)

    # if at least one ID could be converted
    can_be_converted = True if mapping_dict else False
    return dataset_metadata, can_be_converted


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FORMATTING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def format_species(species: str) -> str:
    return "_".join(species.lower().split(" ")[:2])


def format_platform_name(platform_name: str) -> str:
    return platform_name.replace("_", "").replace("-", "").lower()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# METADATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
        "platform_accessions": platform_accessions,
        "summary": dataset_metadata["summary"],
        "title": dataset_metadata["title"],
        "overall_design": additional_metadata["Series"]["Overall-Design"],
        "experiment_types": dataset_metadata["gdsType"],
        "sample_characteristics": list(set(sample_characteristics)),
        "sample_library_strategies": list(set(sample_library_strategies)),
        "sample_library_sources": list(set(sample_library_sources)),
        "sample_descriptions": list(set(sample_descriptions)),
        "sample_titles": list(set(sample_titles)),
        "sample_molecule_types": list(set(sample_molecule_types)),
    }


def parse_metadata(dataset_metadata: dict) -> dict | None:
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
    datasets: list[dict], excluded_accessions_file: Path
) -> list[dict]:
    # parsing list of unwanted accessions
    with open(excluded_accessions_file) as fin:
        excluded_accessions = fin.read().splitlines()
    datasets_to_keep = []
    for dataset in datasets:
        if dataset["Accession"] not in excluded_accessions:
            datasets_to_keep.append(dataset)
    return datasets_to_keep


def species_is_ok(dataset: dict, species: str) -> bool:
    accession = dataset["Accession"]
    # we want datasets only specific to the species we are interested in
    parsed_species_list = dataset["taxon"].split("; ")
    if not parsed_species_list:
        logger.warning(f"Accession {accession} rejected: Could not detect species.")
        return False
    # trying to find our species in the list of species parsed
    for parsed_species in parsed_species_list:
        if format_species(parsed_species) == format_species(species):
            if len(parsed_species_list) > 1:
                logger.info(
                    f"Accession {accession}: multiple species detected = {parsed_species_list}"
                )
            return True
    logger.warning(
        f"Accession {accession} rejected: Found wrong species = {parsed_species_list}"
    )
    return False


def contains_only_rna(molecules_types: list, accession: str) -> bool:
    # we want only GEO series that contain only RNA molecules
    # for other series, they should be superseries contained other series that are being parsed too
    # so anyway, this would lead in duplicates
    if all(["rna" in molecule_type.lower() for molecule_type in molecules_types]):
        return True
    logger.info(f"Accession {accession} rejected: Molecule type(s) = {molecules_types}")
    return False


def contains_proper_experiment_type(
    experiment_types: list, accession: str, platform: str
) -> bool:
    for experiment_type in experiment_types:
        # if at least one experiment type is ok, we keep this dataset
        if GEO_EXPERIMENT_TYPE_TO_PLATFORM.get(experiment_type) == platform:
            return True
    logger.info(
        f"Accession {accession} rejected: Experiment type(s) = {experiment_types}"
    )
    return False


def contains_transcriptomic_source(library_sources: list, accession: str) -> bool:
    # if we have no data about library sources, we just cannot infer
    if not library_sources:
        return True
    # TODO: see how to process series with multiple library sources
    if len(library_sources) > 1:
        return False
    if library_sources[0] in ALLOWED_LIBRARY_SOURCES:
        return True
    logger.warning(f"Accession {accession} rejected: Source(s) = {library_sources}")
    return False


def dataset_is_valid(metadata: dict, platform: str) -> bool:
    accession = metadata["accession"]
    # checking platform
    if not contains_proper_experiment_type(
        metadata["experiment_types"], accession, platform
    ):
        return False

    # checking that library sources fit
    if not contains_transcriptomic_source(
        metadata["sample_library_sources"], accession
    ):
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
        logger.info(
            f"Found keywords: {found_keywords} in accession {metadata['accession']}"
        )
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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING GEO DATASETS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Getting datasets corresponding to species {args.species}")
    dataset_metadata_list = fetch_geo_datasets_for_species(args.species)
    logger.info(
        f"Found {len(dataset_metadata_list)} datasets for species {args.species}"
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FOR DEV PURPOSES / TESTING: RESTRICT TO SPECIFIC ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.accessions:
        logger.info(f"Keeping only accessions {args.accessions}")
        dev_accessions = args.accessions.split(",")
        dataset_metadata_list = [
            d for d in dataset_metadata_list if d["Accession"] in dev_accessions
        ]
        logger.info(
            f"Kept {len(dataset_metadata_list)} datasets for dev / testing purposes"
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXCLUDING UNWANTED ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.excluded_accessions_file:
        logger.info("Excluding unwanted datasets")
        dataset_metadata_list = exclude_unwanted_accessions(
            dataset_metadata_list, args.excluded_accessions_file
        )
        logger.info(
            f"{len(dataset_metadata_list)} datasets remaining after excluding unwanted accessions"
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXCLUDING DATASETS WITH THE WRONG SPECIES
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info("Excluding wrong species")
    tmp_lst = [
        dataset
        for dataset in dataset_metadata_list
        if species_is_ok(dataset, args.species)
    ]

    # checking if all datasets were ok
    if len(tmp_lst) < len(dataset_metadata_list):
        logger.warning(
            f"{len(dataset_metadata_list) - len(tmp_lst)} dataset(s) did not have the correct species!"
        )
        selected_metadata_list = []
    else:
        logger.info("All datasets had the correct species")

    dataset_metadata_list = tmp_lst

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info("Parsing metadata")
    metadata_list = []
    with (
        Pool(processes=args.nb_cpus) as p,
        tqdm(total=len(dataset_metadata_list)) as pbar,
    ):
        for result in p.imap_unordered(parse_metadata, dataset_metadata_list):
            pbar.update()
            pbar.refresh()
            if result is None:
                continue
            metadata_list.append(result)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # CHECKING MOLECULE TYPE / PLATFORM TECHNOLOGIES
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info("Validating datasets")
    filtered_metadata_list = [
        metadata
        for metadata in metadata_list
        if dataset_is_valid(metadata, args.platform)
    ]
    logger.info(
        f"{len(filtered_metadata_list)} datasets remaining after checking technology platform and molecule type"
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FILTERING WITH KEYWORDS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.keywords:
        logger.info(f"Filtering experiments with keywords {args.keywords}")
        func = partial(filter_metadata_with_keywords, keywords=args.keywords)

        selected_metadata_list = []
        with (
            Pool(processes=args.nb_cpus) as p,
            tqdm(total=len(filtered_metadata_list)) as pbar,
        ):
            for result in p.imap_unordered(func, filtered_metadata_list):
                pbar.update()
                pbar.refresh()
                if result is None:
                    continue
                selected_metadata_list.append(result)

        logger.info(
            f"{len(selected_metadata_list)} datasets remaining after filtering with keywords"
        )

    else:
        selected_metadata_list = filtered_metadata_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GETTING METADATA OF SEQUENCING PLATFORMS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info("Getting platform metadata")
    tmp_lst = []
    for selected_metadata_chunk_list in tqdm(
        chunk_list(selected_metadata_list, ENTREZ_CHUNKSIZE)
    ):
        tmp_lst += get_platform_metadata(selected_metadata_chunk_list)

    # checking if platform metadata was found for all datasets
    if len(tmp_lst) < len(selected_metadata_list):
        logger.warning(
            f"Platform metadata could not be retrieved for {len(selected_metadata_list) - len(tmp_lst)} dataset(s)!"
        )
        selected_metadata_list = []
    else:
        logger.info("Platform metadata found for all datasets!")

    # augmenting selected_metadata_list with platform metadata
    selected_metadata_list = tmp_lst

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FILTERING OUT DATASETS FOR WHICH ID MAPPING DOES NOT WORK
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # this cannot be done in parallel because it requires HTTP requests
    logger.info("Checking gene ID mapping issues")
    func = partial(probe_ids_can_be_converted, species=args.species)
    final_metadata_list = []

    with (
        Pool(processes=args.nb_cpus) as p,
        tqdm(total=len(filtered_metadata_list)) as pbar,
    ):
        for metadata, can_be_converted in p.imap_unordered(
            func, selected_metadata_list
        ):
            pbar.update()
            pbar.refresh()
            if can_be_converted:
                final_metadata_list.append(metadata)

    logger.info(
        f"{len(final_metadata_list)} datasets remaining after checking gene ID mapping issues"
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GETTING ACCESSIONS TO DOWNLOAD
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if final_metadata_list:
        logger.info(f"Kept {len(final_metadata_list)} datasets")
        # getting accessions of selected experiments
        selected_accessions = [
            metadata["accession"] for metadata in final_metadata_list
        ]

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
        logger.info(
            f"Writing metadata of filtered datasets to {FILTERED_DATASETS_METADATA_OUTFILE_NAME}"
        )
        df = pd.DataFrame.from_dict(filtered_metadata_list)
        df.to_csv(
            FILTERED_DATASETS_METADATA_OUTFILE_NAME,
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

    rejected_metadata_list = [
        metadata for metadata in metadata_list if metadata not in filtered_metadata_list
    ]
    if rejected_metadata_list:
        logger.info(
            f"Writing metadata of rejected datasets to {REJECTED_DATASETS_METADATA_OUTFILE_NAME}"
        )
        df = pd.DataFrame.from_dict(rejected_metadata_list)
        df.to_csv(
            REJECTED_DATASETS_METADATA_OUTFILE_NAME,
            sep="\t",
            index=False,
            header=True,
        )

    if selected_metadata_list:
        logger.info(
            f"Writing metadata of selected datasets to {SELECTED_DATASETS_METADATA_OUTFILE_NAME}"
        )
        df = pd.DataFrame.from_dict(selected_metadata_list)
        df.to_csv(
            SELECTED_DATASETS_METADATA_OUTFILE_NAME,
            sep="\t",
            index=False,
            header=True,
        )

    if final_metadata_list:
        logger.info(
            f"Writing metadata of selected datasets to {FINAL_DATASETS_METADATA_OUTFILE_NAME}"
        )
        df = pd.DataFrame.from_dict(final_metadata_list)
        df.to_csv(
            FINAL_DATASETS_METADATA_OUTFILE_NAME,
            sep="\t",
            index=False,
            header=True,
        )


if __name__ == "__main__":
    main()
