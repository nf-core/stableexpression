#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from datetime import datetime
from pathlib import Path
from urllib.request import urlretrieve

import pandas as pd
import requests
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
SPECIES_INFO_EXT = "info/genomes/taxonomy/{species}"
ENSEMBL_API_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}
STOP_RETRY_AFTER_DELAY = 600

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
# REQUESTS
##################################################################
##################################################################


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def parse_page_data(url: str) -> BeautifulSoup:
    page = requests.get(url)
    page.raise_for_status()
    return BeautifulSoup(page.content, "html.parser")


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(taxid: str | int):
    taxons = [str(taxid)]
    data = {"taxons": taxons}
    response = requests.post(NCBI_TAXONOMY_API_URL, headers=NCBI_API_HEADERS, json=data)
    response.raise_for_status()
    return response.json()


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_get_request_to_ensembl(url: str) -> list[dict]:
    logger.info(f"Sending GET request to {url}")
    response = requests.get(url, headers=ENSEMBL_API_HEADERS)
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
    result = send_request_to_ncbi_taxonomy(species)
    if len(result["taxonomy_nodes"]) > 1:
        raise ValueError(f"Multiple taxids for species {species}")
    metadata = result["taxonomy_nodes"][0]
    if "taxonomy" not in metadata:
        raise ValueError(f"Could not find taxonomy results for species {species}")
    return int(metadata["taxonomy"]["tax_id"])


def get_species_division(species_taxid: int) -> str:
    url = ENSEMBL_REST_SERVER + SPECIES_INFO_EXT.format(species=str(species_taxid))
    data = send_get_request_to_ensembl(url)
    if len(data) == 0:
        raise ValueError(f"No division found for species Taxon ID {species_taxid}")
    elif len(data) > 1:
        logger.warning(
            f"Multiple divisions found for species Taxon ID {species_taxid}. Keeping the first one."
        )
    return data[0]["division"]


def get_species_category(species: str) -> str:
    ncbi_formated_species_name = format_species_name_for_ncbi_taxonomy(species)
    species_taxid = get_species_taxid(ncbi_formated_species_name)
    division = get_species_division(species_taxid)
    return ENSEMBL_DIVISION_TO_FOLDER[division]


def get_division_url(species: str) -> str:
    category = get_species_category(species)
    if category == "vertebrates":
        return ENSEMBL_VERTEBRATES_BASE_URL
    else:
        return ENSEMBL_GENOMES_BASE_URL.format(category)


def format_species_name_for_ensembl(species: str) -> str:
    return species.replace(" ", "_").lower()


def format_species_name_for_ncbi_taxonomy(species: str) -> str:
    return species.replace("_", " ").lower()


def parse_last_modified_date(dt_string: str) -> datetime | None:
    try:
        return datetime.strptime(dt_string, "%Y-%m-%d %H:%M")
    except ValueError:
        return None


def get_candidate_species_folders(
    species: str, url: str, first_level: bool = True
) -> list[dict]:
    soup = parse_page_data(url)
    species_url_records = []

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
            folder_url = f"{url}{folder.text}"
            if folder.text.startswith(species):
                d = {
                    "date": last_modified_date,
                    "url": folder_url,
                    "name": folder.text.rstrip("/"),
                }
                species_url_records.append(d)
                print(folder.text)
            elif folder.text.endswith("_collection/"):
                species_url_records += get_candidate_species_folders(
                    species, folder_url, first_level=False
                )
            else:
                continue

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

    species = format_species_name_for_ensembl(args.species)
    division_url = get_division_url(species)
    logger.info(f"Searching for the right folder in {division_url}")

    species_url_records = get_candidate_species_folders(species, division_url)
    if not species_url_records:
        raise ValueError(f"No species folder found for {species}")

    annotation_folder_url = get_current_annotation_folder(species_url_records, species)
    logger.info(f"Found current annotation folder: {annotation_folder_url}")

    annotation_file = get_annotation_file(annotation_folder_url)

    annotation_full_url = annotation_folder_url + annotation_file
    logger.info(f"Found annotation URL: {annotation_full_url}.\nDownloading...")

    download_file(annotation_full_url, annotation_file)
    logger.info("Done")


if __name__ == "__main__":
    main()
