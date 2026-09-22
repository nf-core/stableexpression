#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
import httpx
import pandas as pd
from typing import ClassVar
from pathlib import Path

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
from ncbi_annotation_manager import NCBIAnnotationManager

logger = logging.getLogger(__name__)
logging.getLogger("httpx").setLevel(logging.ERROR)


STOP_RETRY_AFTER_DELAY = 120

ENSEMBL_API_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}


##################################################################
##################################################################
# REQUEST FUNCTIONS
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


##################################################################
##################################################################
# PARSING FUNCTIONS
##################################################################
##################################################################

def parse_last_modified_date(dt_string: str) -> datetime | None:
    try:
        return datetime.strptime(dt_string, "%Y-%m-%d %H:%M")
    except ValueError:
        return None


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


##################################################################
##################################################################
# CLASS
##################################################################
##################################################################


class EnsemblAnnotationManager:

    ENSEMBL_REST_SERVER: ClassVar[str] = "https://rest.ensembl.org/"
    SPECIES_INFO_BASE_ENDPOINT: ClassVar[str] = "info/genomes/taxonomy/{species}"
    TAXONOMY_NAME_ENDPOINT: ClassVar[str] = "taxonomy/name/{species}"

    ENSEMBL_DIVISION_TO_FOLDER: ClassVar[dict[str, str]] = {
        "EnsemblPlants": "plants",
        "EnsemblVertebrates": "vertebrates",
        "EnsemblMetazoa": "metazoa",
        "EnsemblFungi": "fungi",
        "EnsemblBacteria": "bacteria",
        "EnsemblProtists": "protists",
    }

    ENSEMBL_GENOMES_BASE_URL: ClassVar[str] = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/{}/gff3/"
    ENSEMBL_VERTEBRATES_BASE_URL: ClassVar[str] = "https://ftp.ensembl.org/pub/current/gff3/"

    species: str
    species_taxid: int
    division: str
    assembly_names: list[str]
    division_url: str


    def __init__(self, species: str):
        self.species = species
        self.species_taxid = self.get_species_taxid()
        logger.info(f"[Ensembl] :: Got species taxid: {self.species_taxid}")
        self.division, self.assembly_names = self.get_species_division_and_candidate_folders()
        logger.info(f"[Ensembl] :: Got division: {self.division}")
        self.division_url = self.get_division_url()


    def get_candidate_folders(self):
        return self._get_candidate_species_folders(self.division_url)


    def get_species_taxid(self) -> int:
        try:
            return self.get_species_taxid_from_ensembl()
        except Exception as e:
            logger.error(
                f"Could not get species taxid for species {self.species} using the Ensembl REST API: {e}.\nTrying NCBI taxonomy."
            )
            return NCBIAnnotationManager.get_species_taxid(self.species)


    def get_taxonomy_metadata(self):
        formatted_species = self.species.replace(" ", "_").lower()
        url = self.ENSEMBL_REST_SERVER + self.TAXONOMY_NAME_ENDPOINT.format(species=formatted_species)
        return send_get_request_to_ensembl(url)


    def get_species_taxid_from_ensembl(self) -> int:
        data: list[dict] = self.get_taxonomy_metadata()
        if len(data) == 0:
            raise ValueError(f"No species found for species {self.species}")
        elif len(data) > 1:
            logger.warning(
                f"Multiple species found for species {self.species}. Keeping the first one."
            )
        species_data = data[0]
        if "id" not in species_data:
            raise ValueError(
                f"Could not find taxid for species {self.species}. Data collected: {species_data}"
            )
        return species_data["id"]


    def get_species_info(self):
        url = self.ENSEMBL_REST_SERVER + self.SPECIES_INFO_BASE_ENDPOINT.format(
            species=str(self.species_taxid)
        )
        return send_get_request_to_ensembl(url)


    def get_species_division_and_candidate_folders(self) -> tuple[str, list[str]]:
        data: list[dict] = self.get_species_info()
        if len(data) == 0:
            raise ValueError(f"No division found for species Taxon ID {self.species_taxid}")
        found_divisions = list({d["division"] for d in data})
        # this should not happen (and if it does, it's an issue on Ensembl's side)
        if not found_divisions:
            raise ValueError(f"Could not find any division for species {self.species_taxid}...")
        # we should never have multiple possible divisions for a single species
        # it is like if a species belonged to multiple kingdoms at the same time...
        if len(found_divisions) > 1:
            raise ValueError(f"Multiple divisions found for species Taxon ID {self.species_taxid}: {found_divisions}.")
        # there should be only one division
        found_division = found_divisions[0]
        # taking all assembly names
        assembly_names = list({d.get("name", "") for d in data})
        return found_division, assembly_names


    def get_division_url(self) -> str:
        """
        In Ensembl, species are separated into divisions.
        Returns the URL for the division of the given species.
        """
        if self.division == "vertebrates":
            return self.ENSEMBL_VERTEBRATES_BASE_URL
        else:
            division_folder = self.ENSEMBL_DIVISION_TO_FOLDER[self.division]
            return self.ENSEMBL_GENOMES_BASE_URL.format(division_folder)


    def _get_candidate_species_folders(self, url: str, first_level: bool = True) -> list[str]:
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
                if folder_name.startswith(self.species) or folder_name in self.assembly_names:
                    folder_urls.append(folder_url)
                elif folder_name.endswith("_collection/"):
                    collection_folder_urls += self._get_candidate_species_folders(folder_url, first_level=False)
                else:
                    continue

        # if no first-level folder found
        if not folder_urls:
            # if collection folders were found at >= second level, taking those ones as fallback
            if collection_folder_urls:
                logger.warning(f"No first-level folder found for {self.species} at {url}. Taking collection folders as fallback.")
                return collection_folder_urls
            else:
                return []

        # if first-level folders were found, keeping only those ones
        return folder_urls


    @staticmethod
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


    @staticmethod
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
