#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
from tqdm import tqdm
from multiprocessing import Pool
from Bio import Entrez
from pathlib import Path

# from random import sample
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
import logging
from requests.exceptions import HTTPError, ConnectionError

from natural_language_utils import keywords_in_fields
from gprofiler_utils import chunk_list

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# set a custom writable directory before any Entrez operations
# mandatory for running the script in an apptainer container
# Entrez.Parser.Parser.directory("/tmp/biopython")

ACCESSION_OUTFILE_NAME = "accessions.tsv"
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

# NCBI_API_BASE_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc={accession}"
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


def augment_with_platform_metadata(
    datasets: list[dict],
) -> tuple[list[dict], list[dict]]:
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
    # one single request to NCBI for all platform accessions
    # we extract the platform accessions to allow better parsing afterwards
    acc_to_metadata = {
        platform_metadata["Accession"]: platform_metadata
        for platform_metadata in fetch_geo_platform_data(platform_accessions)
    }

    # adding the platform metadata to the corresponding metadata
    issues = []
    augmented_metadata_list = []
    for dataset in datasets:
        accession = dataset["accession"]
        platform_accessions = dataset["platform_accessions"]
        dataset["platform_metadata"] = []

        if not platform_accessions:
            issues.append({"accession": accession, "reason": "NO PLATFORM ACCESSIONS"})
            continue

        for platform_accession in platform_accessions:
            # filtering out cases where the platform metadata is not available
            if platform_accession not in acc_to_metadata:
                continue
            # augmenting metadata with platform metadata
            dataset["platform_metadata"].append(acc_to_metadata[platform_accession])

        # getting list of platform taxon
        platforms_taxons = [
            platform_metadata.get("taxon")
            for platform_metadata in dataset["platform_metadata"]
            if platform_metadata.get("taxon") is not None
        ]

        # checking if there is one single platform taxon
        # otherwise, checking the dataset
        if not platforms_taxons:
            logger.warning(f"No taxon found for dataset {accession}")
            issues.append({"accession": accession, "reason": "NO PLATFORM TAXON"})
            continue
        elif len(platforms_taxons) > 1:
            logger.warning(
                f"Multiple taxons for dataset {accession}: {platforms_taxons}"
            )
            issues.append(
                {
                    "accession": accession,
                    "reason": f"MULTIPLE PLATFORM TAXONS: {platforms_taxons}",
                }
            )
            continue

        dataset["platform_taxon"] = platforms_taxons[0]
        augmented_metadata_list.append(dataset)

    return augmented_metadata_list, issues


"""
def download_platform_datatable(ftp_link: str, platform_accession: str) -> Path | None:
    filename = f"soft/{platform_accession}_family.soft.gz"
    ftp_url = ftp_link + filename
    output_file = Path(PLATFORM_SOFT_TMPDIR) / f"{platform_accession}.gz"
    download_file_at_url(ftp_url, output_file)
    return output_file
"""

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
"""

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
        "taxon": dataset_metadata["taxon"],
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


def check_species_issues(parsed_species_list: list, species: str) -> dict:
    # trying to find our species in the list of species parsed
    for parsed_species in parsed_species_list:
        if format_species(parsed_species) == format_species(species):
            return {}
    return {"parsed_species": parsed_species_list}


def check_molecule_type_issues(molecules_types: list) -> dict:
    # we want only GEO series that contain only RNA molecules
    # for other series, they should be superseries contained other series that are being parsed too
    # so anyway, this would lead in duplicates
    if all(["rna" in molecule_type.lower() for molecule_type in molecules_types]):
        return {}
    return {"molecule_types": molecules_types}


def check_experiment_type_issues(experiment_types: list | str, platform: str) -> dict:
    experiment_types = (
        experiment_types if isinstance(experiment_types, list) else [experiment_types]
    )
    for experiment_type in experiment_types:
        # if at least one experiment type is ok, we keep this dataset
        if GEO_EXPERIMENT_TYPE_TO_PLATFORM.get(experiment_type) == platform:
            return {}
    return {"experiment_types": experiment_types}


def check_source_issues(library_sources: list) -> dict:
    # if we have no data about library sources, we just cannot infer
    if not library_sources:
        return {}
    if len(library_sources) == 1 and library_sources[0] in ALLOWED_LIBRARY_SOURCES:
        return {}
    # TODO: see how to process series with multiple library sources
    return {"library_sources": library_sources}


def search_keywords(dataset: dict, keywords: list[str]) -> tuple[list, dict]:
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
        return found_keywords, {}
    else:
        return [], {"accession": accession, "keywords_found": False}


def check_dataset(
    dataset: dict, species: str, platform: str | None, keywords: list[str] | None
) -> tuple[list, dict]:
    accession = dataset["accession"]
    parsed_species_list = dataset["taxon"].split("; ")
    experiment_types = dataset["experiment_types"]
    library_sources = dataset["sample_library_sources"]
    molecules_types = dataset["sample_molecule_types"]

    # checking species
    issue_dict = check_species_issues(parsed_species_list, species)

    # checking platform
    if platform is not None:
        platform_issue_dict = check_experiment_type_issues(experiment_types, platform)
        issue_dict |= platform_issue_dict

    # checking that library sources fit
    transcriptomic_issue_dict = check_source_issues(library_sources)
    issue_dict |= transcriptomic_issue_dict

    # checking that all molecule types are RNA
    moltype_issue_dict = check_molecule_type_issues(molecules_types)
    issue_dict |= moltype_issue_dict

    found_keywords = []
    if keywords:
        found_keywords, keyword_issue_dict = search_keywords(dataset, keywords)
        issue_dict |= keyword_issue_dict

    if issue_dict:
        rejection_dict = {"accession": accession, "reasons": issue_dict}
    else:
        rejection_dict = {}

    return found_keywords, rejection_dict


def export_dataset_metadatas(
    datasets: list[dict], output_file: str, clean_columns: bool = True
):
    if datasets:
        df = pd.DataFrame.from_dict(datasets)
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
    # EXCLUDING UNWANTED ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.excluded_accessions_file:
        logger.info("Excluding unwanted datasets")
        datasets = exclude_unwanted_accessions(datasets, args.excluded_accessions_file)
        logger.info(
            f"{len(datasets)} datasets remaining after excluding unwanted accessions"
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Parsing metadata for {len(datasets)} datasets")
    augmented_datasets = []
    with (
        Pool(processes=args.nb_cpus) as p,
        tqdm(total=len(datasets)) as pbar,
    ):
        for result in p.imap_unordered(parse_metadata, datasets):
            pbar.update()
            pbar.refresh()
            if result is None:
                continue
            augmented_datasets.append(result)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # CHECKING DATASET METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Validating {len(augmented_datasets)} datasets")
    selected_datasets = []
    rejected_datasets = []
    for dataset in tqdm(augmented_datasets):
        found_keywords, rejection_dict = check_dataset(
            dataset, args.species, args.platform, args.keywords
        )
        if rejection_dict:
            rejected_datasets.append(rejection_dict)
        else:
            if found_keywords:
                dataset["found_keywords"] = found_keywords
            selected_datasets.append(dataset)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GETTING METADATA OF SEQUENCING PLATFORMS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Getting platform metadata for {len(selected_datasets)} datasets")
    selected_datasets_chunks = chunk_list(
        selected_datasets, PLATFORM_METADATA_CHUNKSIZE
    )
    # resetting selecting datasets
    selected_datasets = []
    for selected_datasets_chunk in tqdm(selected_datasets_chunks):
        augmented_datasets, issues = augment_with_platform_metadata(
            selected_datasets_chunk
        )
        selected_datasets += augmented_datasets
        rejected_datasets += issues

    """
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FILTERING OUT DATASETS FOR WHICH ID MAPPING DOES NOT WORK
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # this cannot be done in parallel because it requires HTTP requests
    logger.info(
        f"Checking gene ID mapping issues for {len(platform_augmented_dataset_metadata_list)} datasets"
    )
    func = partial(probe_ids_can_be_converted, species=args.species)
    final_metadata_list = []

    with (
        Pool(processes=args.nb_cpus) as p,
        tqdm(total=len(platform_augmented_dataset_metadata_list)) as pbar,
    ):
        for metadata, can_be_converted in p.imap_unordered(
            func, platform_augmented_dataset_metadata_list
        ):
            pbar.update()
            pbar.refresh()
            if can_be_converted:
                final_metadata_list.append(metadata)

    export_filtered_out_datasets_if_any(
        platform_augmented_dataset_metadata_list,
        final_metadata_list,
        GENE_ID_MAPPING_ISSUES_DATASETS_METADATA_OUTFILE_NAME,
        "gene id mapping",
    )
    """

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXPORTING ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info(f"Kept {len(selected_datasets)} datasets")
    # getting accessions of selected experiments
    selected_accessions = [
        {"accession": dataset["accession"], "platform_taxon": dataset["platform_taxon"]}
        for dataset in selected_datasets
    ]
    export_dataset_metadatas(selected_accessions, ACCESSION_OUTFILE_NAME)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXPORTING DATASETS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    export_dataset_metadatas(augmented_datasets, SPECIES_DATASETS_OUTFILE_NAME)
    export_dataset_metadatas(selected_datasets, SELECTED_DATASETS_OUTFILE_NAME)
    export_dataset_metadatas(
        rejected_datasets, REJECTED_DATASETS_OUTFILE_NAME, clean_columns=False
    )


if __name__ == "__main__":
    main()
