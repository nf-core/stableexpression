#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from datetime import datetime
from urllib.request import urlretrieve

import httpx
import pandas as pd
from bs4 import BeautifulSoup
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

GENE_IDS_CHUNKSIZE = 50  # max allowed by Ensembl REST API

ENSEMBL_REST_SERVER = "https://rest.ensembl.org/"
SPECIES_INFO_BASE_ENDPOINT = "info/genomes/taxonomy/{species}"
TAXONOMY_NAME_ENDPOINT = "taxonomy/name/{species}"
ENSEMBL_API_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}
STOP_RETRY_AFTER_DELAY = 120

NCBI_TAXONOMY_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}

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
        dest="species",
        required=True,
        help="Species name",
    )
    return parser.parse_args()


##################################################################
##################################################################
# httpx
##################################################################
##################################################################


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def parse_page_data(url: str) -> BeautifulSoup:
    page = httpx.get(url)
    page.raise_for_status()
    return BeautifulSoup(page.content, "html.parser")


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(taxid: str | int):
    """
    Sends a POST request to the NCBI taxonomy API to retrieve taxonomic information for the given taxid.
    """
    logger.info(f"Sending POST request to {NCBI_TAXONOMY_API_URL}")
    taxons = [str(taxid)]
    data = {"taxons": taxons}
    response = httpx.post(NCBI_TAXONOMY_API_URL, headers=NCBI_API_HEADERS, json=data)
    response.raise_for_status()
    return response.json()


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_get_request_to_ensembl(url: str) -> list[dict]:
    """
    Sends a GET request to the Ensembl API to retrieve data from the given URL.
    """
    logger.info(f"Sending GET request to {url}")
    response = httpx.get(url, headers=ENSEMBL_API_HEADERS)
    if response.status_code == 200:
        response.raise_for_status()
    else:
        raise RuntimeError(
            f"Failed to retrieve data: encountered error {response.status_code}"
        )
    return response.json()


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def download_file(url: str, output_path: str):
    try:
        urlretrieve(url, output_path)
    except Exception as e:
        logger.error(f"Failed to download file from {url}: {e}")
        raise


##################################################################
##################################################################
# PARSING
##################################################################
##################################################################


def get_species_taxid(species: str) -> int:
    try:
        return get_species_taxid_from_ensembl(species)
    except Exception as e:
        logger.error(
            f"Could not get species taxid for species {species} using the Ensembl REST API: {e}.\nTrying NCBI taxonomy."
        )
        ncbi_formated_species_name = format_species_name_for_ncbi_taxonomy(species)
        return get_species_taxid_from_ncbi(ncbi_formated_species_name)


def format_species_name_for_ensembl(species: str) -> str:
    return species.replace(" ", "_").lower()


def format_species_name_for_ncbi_taxonomy(species: str) -> str:
    return species.replace("_", " ").lower()


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


def get_species_taxid_from_ncbi(species: str) -> int:
    formatted_species = format_species_name_for_ncbi_taxonomy(species)
    result = send_request_to_ncbi_taxonomy(formatted_species)
    if len(result["taxonomy_nodes"]) > 1:
        raise ValueError(f"Multiple taxids for species {species}")
    metadata = result["taxonomy_nodes"][0]
    if "taxonomy" not in metadata:
        raise ValueError(f"Could not find taxonomy results for species {species}")
    return int(metadata["taxonomy"]["tax_id"])


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
        logger.error(f"Could not find any division for species {species_taxid}...")
        sys.exit(100)
    # we should never have multiple possible divisions for a single species
    # it is like if a species belonged to multiple kingdoms at the same time...
    if len(found_divisions) > 1:
        logger.error(f"Multiple divisions found for species Taxon ID {species_taxid}: {found_divisions}.")
        sys.exit(100)
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
) -> list[dict]:
    """
    Get the content of the url corresponding to the species division and parse it using BeautifulSoup.
    Get all folders
    """
    soup = parse_page_data(url)
    species_url_records = []
    species_url_collection_records = []

    # adding progress bar only at the first level
    iterator = tqdm(soup.find_all("tr")) if first_level else soup.find_all("tr")
    for item in iterator:
        # all line sections
        line_sections = list(item.find_all("td"))
        # all folders of interest have an associated date
        if len(line_sections) < 2:
            continue

        folder_name_section = line_sections[1]
        date_section = line_sections[2]
        last_modified_date = parse_last_modified_date(date_section.text.strip())

        for folder in folder_name_section.find_all("a"):
            folder_name = folder.text
            folder_url = f"{url}{folder_name}"
            # getting all folders that either
            # start with the species name
            # are in the list of assembly names
            if folder_name.startswith(species) or folder_name in assembly_names:
                d = {
                    "date": last_modified_date,
                    "url": folder_url,
                    "name": folder_name.rstrip("/"),
                }
                species_url_records.append(d)
            elif folder_name.endswith("_collection/"):
                species_url_collection_records += get_candidate_species_folders(
                    species, assembly_names, folder_url, first_level=False
                )
            else:
                continue

    # if no first-level folder found
    if not species_url_records:
        # if collection folders were found at >= second level, taking those ones as fallback
        if species_url_collection_records:
            logger.warning(f"No first-level folder found for {species} at {url}. Taking collection folders as fallback.")
            return species_url_collection_records
        else:
            raise ValueError(f"No species folder found for {species} at {url}")

    # if first-level folders were found, keeping only those ones
    return species_url_records


def get_main_folder_url(records: list[dict], species: str) -> str | None:
    main_folder_url = None
    for record in records:
        if record["name"] == species:
            main_folder_url = record["url"]
            break
    return main_folder_url


def get_last_modified_folder_url(records: list[dict]) -> str:
    df = pd.DataFrame.from_dict(records)
    df.sort_values(by="date", ascending=False, inplace=True)
    return df.iloc[0]["url"]


def get_current_annotation_folder(records: list[dict], species: str) -> str:
    main_folder_url = get_main_folder_url(records, species)
    if main_folder_url is not None:
        return main_folder_url

    logger.info(
        "Could not find a folder having the species as name. Checking for gca folders."
    )
    gca_records = [
        record for record in records if record["name"].startswith(f"{species}_gca")
    ]
    if gca_records:
        return get_last_modified_folder_url(gca_records)

    logger.info(
        "Could not find a folder having the species as name. Getting the last modified one."
    )
    return get_last_modified_folder_url(records)


def parse_size(size_str):
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
    multipliers = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4, "P": 1024**5}

    # Check if last character is a unit
    if size_str[-1] in multipliers:
        number = float(size_str[:-1])
        multiplier = multipliers[size_str[-1]]
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

    df = pd.DataFrame.from_dict(file_records)

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


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    args = parse_args()

    species_taxid = get_species_taxid(args.species)
    logger.info(f"Got species taxid: {species_taxid}")

    division, assembly_names = get_species_division_and_candidate_folders(species_taxid)
    logger.info(f"Got division: {division}")

    logger.info(f"Fetching division name for {args.species}")
    division_url = get_division_url(division)

    logger.info(f"Searching for the right folder in {division_url}")
    species_url_records = get_candidate_species_folders(args.species, assembly_names, division_url)
    if not species_url_records:
        raise ValueError(f"No species folder found for {args.species}")

    annotation_folder_url = get_current_annotation_folder(species_url_records, args.species)
    logger.info(f"Found current annotation folder: {annotation_folder_url}")

    annotation_file = get_annotation_file(annotation_folder_url)

    annotation_full_url = annotation_folder_url + annotation_file
    logger.info(f"Found annotation URL: {annotation_full_url}.\nDownloading...")

    download_file(annotation_full_url, annotation_file)
    logger.info("Done")


if __name__ == "__main__":
    main()
