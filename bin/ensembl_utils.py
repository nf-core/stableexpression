#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
import httpx
import pandas as pd

from datetime import datetime
from tqdm import tqdm
from urllib.request import urlretrieve

from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)
from bs4 import BeautifulSoup

import ncbi_datasets_utils

logger = logging.getLogger(__name__)
logging.getLogger("httpx").setLevel(logging.ERROR)


ENSEMBL_REST_SERVER = "https://rest.ensembl.org/"
SPECIES_INFO_BASE_ENDPOINT = "info/genomes/taxonomy/{species}"
TAXONOMY_NAME_ENDPOINT = "taxonomy/name/{species}"
ENSEMBL_API_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}
ENSEMBL_DIVISION_TO_FOLDER = {
    "EnsemblPlants": "plants",
    "EnsemblVertebrates": "vertebrates",
    "EnsemblMetazoa": "metazoa",
    "EnsemblFungi": "fungi",
    "EnsemblBacteria": "bacteria",
    "EnsemblProtists": "protists",
}

ENSEMBL_GENOMES_BASE_URL = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/{}/gff3/"
ENSEMBL_VERTEBRATES_BASE_URL = "https://ftp.ensembl.org/pub/current/gff3/"

STOP_RETRY_AFTER_DELAY = 120


##################################################################
##################################################################
# FUNCTIONS
##################################################################
##################################################################


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_get_request_to_ensembl(url: str) -> list[dict]:
    """
    Sends a GET request to the Ensembl API to retrieve data from the given URL.
    """
    with httpx.Client() as client:
        response = client.get(url, headers=ENSEMBL_API_HEADERS)
        if response.status_code == 200:
            response.raise_for_status()
        else:
            raise RuntimeError(
                f"Failed to retrieve data: encountered error {response.status_code}"
            )
        return response.json()


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def parse_page_data(url: str) -> BeautifulSoup:
    with httpx.Client() as client:
        page = client.get(url)
        page.raise_for_status()
        return BeautifulSoup(page.content, "html.parser")


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def download_file(url: str, output_path: Path):
    try:
        urlretrieve(url, str(output_path))
    except Exception as e:
        logger.error(f"Failed to download file from {url}: {e}")
        raise


def get_species_taxid(species: str) -> int:
    try:
        return get_species_taxid_from_ensembl(species)
    except Exception as e:
        logger.error(
            f"Could not get species taxid for species {species} using the Ensembl REST API: {e}.\nTrying NCBI taxonomy."
        )
        ncbi_formated_species_name = ncbi_datasets_utils.format_species_name(species)
        return ncbi_datasets_utils.get_species_taxid(ncbi_formated_species_name)


def format_species_name_for_ensembl(species: str) -> str:
    return species.replace(" ", "_").lower()


def get_species_taxid_from_ensembl(species: str) -> int:
    formatted_species = format_species_name_for_ensembl(species)
    url = ENSEMBL_REST_SERVER + TAXONOMY_NAME_ENDPOINT.format(species=formatted_species)
    data = send_get_request_to_ensembl(url)
    if len(data) == 0:
        raise ValueError(f"No species found for species {species}")
    elif len(data) > 1:
        logger.warning(
            f"Multiple species found for species {species}. Keeping the first one."
        )
    species_data = data[0]
    if "id" not in species_data:
        raise ValueError(
            f"Could not find taxid for species {species}. Data collected: {species_data}"
        )
    return species_data["id"]


def get_species_division_and_candidate_folders(species_taxid: int) -> tuple[str, list[str]]:
    url = ENSEMBL_REST_SERVER + SPECIES_INFO_BASE_ENDPOINT.format(
        species=str(species_taxid)
    )
    data: list[dict] = send_get_request_to_ensembl(url)
    if len(data) == 0:
        raise ValueError(f"No division found for species Taxon ID {species_taxid}")
    found_divisions = list({d["division"] for d in data})
    # this should not happen (and if it does, it's an issue on Ensembl's side)
    if not found_divisions:
        raise ValueError(f"Could not find any division for species {species_taxid}...")
    # we should never have multiple possible divisions for a single species
    # it is like if a species belonged to multiple kingdoms at the same time...
    if len(found_divisions) > 1:
        raise ValueError(f"Multiple divisions found for species Taxon ID {species_taxid}: {found_divisions}.")
    # there should be only one division
    found_division = found_divisions[0]
    # taking all assembly names
    assembly_names = list({d.get("name", "") for d in data})
    return found_division, assembly_names


def get_division_url(division: str) -> str:
    """
    In Ensembl, species are separated into divisions.
    Returns the URL for the division of the given species.
    """
    if division == "vertebrates":
        return ENSEMBL_VERTEBRATES_BASE_URL
    else:
        division_folder = ENSEMBL_DIVISION_TO_FOLDER[division]
        return ENSEMBL_GENOMES_BASE_URL.format(division_folder)


def parse_last_modified_date(dt_string: str) -> datetime | None:
    try:
        return datetime.strptime(dt_string, "%Y-%m-%d %H:%M")
    except ValueError:
        return None


def get_candidate_species_folders(
    species: str, assembly_names: list[str], url: str, first_level: bool = True
) -> list[str]:
    """
    Get the content of the url corresponding to the species division and parse it using BeautifulSoup.
    Get all folders
    """
    soup = parse_page_data(url)
    folder_urls = []
    collection_folder_urls = []

    # adding progress bar only at the first level
    iterator = tqdm(soup.find_all("tr")) if first_level else soup.find_all("tr")
    for item in iterator:
        # all line sections
        line_sections = list(item.find_all("td"))
        # all folders of interest have an associated date
        if len(line_sections) < 2:
            continue

        folder_name_section = line_sections[1]
        for folder in folder_name_section.find_all("a"):
            folder_name = folder.text
            folder_url = f"{url}{folder_name}"
            # getting all folders that either
            # start with the species name
            # are in the list of assembly names
            if folder_name.startswith(species) or folder_name in assembly_names:
                folder_urls.append(folder_url)
            elif folder_name.endswith("_collection/"):
                collection_folder_urls += get_candidate_species_folders(
                    species, assembly_names, folder_url, first_level=False
                )
            else:
                continue

    # if no first-level folder found
    if not folder_urls:
        # if collection folders were found at >= second level, taking those ones as fallback
        if collection_folder_urls:
            logger.warning(f"No first-level folder found for {species} at {url}. Taking collection folders as fallback.")
            return collection_folder_urls
        else:
            return []

    # if first-level folders were found, keeping only those ones
    return folder_urls


def parse_size(size_str: str) -> int:
    """
    Convert size strings like '902K', '4.1M', '5G' to bytes.

    Parameters:
    -----------
    size_str : str
        Size string with suffix (K, M, G, T, etc.)

    Returns:
    --------
    int : size in bytes
    """
    size_str = size_str.strip().upper()
    # Define multipliers
    MULTIPLIERS = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4, "P": 1024**5}
    # Check if last character is a unit
    if size_str[-1] in MULTIPLIERS:
        number = float(size_str[:-1])
        multiplier = MULTIPLIERS[size_str[-1]]
        return int(number * multiplier)
    else:
        # No suffix, assume it's already in bytes
        return int(float(size_str))


def get_annotation_file(url: str) -> str:
    soup = parse_page_data(url)
    file_records = []

    for item in soup.find_all("tr"):
        # all line sections
        line_sections = list(item.find_all("td"))
        if len(line_sections) < 4:
            continue

        file = line_sections[1].text.strip()
        if not file.endswith(".gff3.gz"):
            continue

        d = {
            "file": file,
            "date": parse_last_modified_date(line_sections[2].text.strip()),
            "size": parse_size(line_sections[3].text.strip()),
        }
        file_records.append(d)

    if not file_records:
        raise ValueError("No annotation files found")

    df = pd.DataFrame(file_records)

    # keeping the biggest annotation
    max_size_df = df.loc[
        [df["size"].idxmax()]
    ]  # double brackets to keep it as a DataFrame
    if len(max_size_df) == 1:
        return max_size_df["file"].iloc[0]

    # if multiple files with the same size, return the most recent
    most_recent_df = max_size_df.loc[
        [max_size_df["date"].idxmax()]
    ]  # double brackets to keep it as a DataFrame
    if len(most_recent_df) == 1:
        return max_size_df["file"].iloc[0]

    # if still multiple files, return the first one
    # remove the one ending with 'chr.gff3.gz' if it exists
    if max_size_df["file"].str.endswith("chr.gff3.gz").any():
        max_size_df = max_size_df[~max_size_df["file"].str.endswith("chr.gff3.gz")]
    return max_size_df["file"].iloc[0]
