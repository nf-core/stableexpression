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
    retry_if_exception
)
from bs4 import BeautifulSoup
from ncbi_annotation_manager import NCBIAnnotationManager

logger = logging.getLogger(__name__)
logging.getLogger("httpx").setLevel(logging.ERROR)

TIMEOUT = 600

STOP_RETRY_AFTER_DELAY = 120

ENSEMBL_API_HEADERS = {
    "Content-type": "application/json"
}


##################################################################
##################################################################
# REQUEST FUNCTIONS
##################################################################
##################################################################

def is_retryable(exception: BaseException) -> bool:
    """Retry everything except a 404 HTTPStatusError."""
    if isinstance(exception, httpx.HTTPStatusError):
        return exception.response.status_code != 404
    return True

@retry(
    retry=retry_if_exception(is_retryable),
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_get_request_to_ensembl(url: str) -> list[dict]:
    """
    Sends a GET request to the Ensembl API to retrieve data from the given URL.
    """
    with httpx.Client(timeout=TIMEOUT) as client:
        response = client.get(url, headers=ENSEMBL_API_HEADERS)
        if response.status_code == 200:
            response.raise_for_status()
        else:
            raise RuntimeError(
                f"Failed to retrieve data: encountered error {response.status_code}"
            )
        return response.json()


@retry(
    retry=retry_if_exception(is_retryable),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def parse_page_data(url: str) -> BeautifulSoup:
    with httpx.Client(timeout=TIMEOUT) as client:
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
        "EnsemblMetazoa": "metazoa",
        "EnsemblFungi": "fungi",
        "EnsemblBacteria": "bacteria",
        "EnsemblProtists": "protists",
    }
    ENSEMBL_VERTEBRATE_DIVISION: ClassVar[str] = "EnsemblVertebrates"

    ENSEMBL_GENOMES_BASE_URL: ClassVar[str] = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/{}/gff3/"
    ENSEMBL_VERTEBRATES_BASE_URL: ClassVar[str] = "https://ftp.ensembl.org/pub/current/gff3/"
    ENSEMBL_ORGANISMS_BASE_URL: ClassVar[str] = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
    ENSEMBL_ORGANISMS_ANNOTATION_FILENAME: ClassVar[str] = "genes.gff3.gz"
    ENSEMBL_ORGANISMS_EXCLUDED_FOLDERS: ClassVar[list[str]] = [
        "Parent Directory"
    ]

    species: str
    species_taxid: int
    target_folder: Path


    def __init__(self, species: str, local_folder: str):
        self.species = species
        self.target_folder = Path(local_folder)
        self.target_folder.mkdir(parents=True, exist_ok=True)
        self.species_taxid = self.get_species_taxid()
        logger.info(f"[Ensembl] :: Got species taxid: {self.species_taxid}")



    def get_ensembl_genomes_annotations(self) -> list[Path]:
        species_info = self.get_species_info()

        if len(species_info) == 0:
            raise ValueError(f"No division found for species Taxon ID {self.species_taxid}")

        division = self.get_species_division(species_info)
        assembly_names = self.get_assembly_names(species_info)
        logger.info(f"[Ensembl] :: Got division: {division}")

        division_url = self.get_division_url(division)
        candidate_folder_urls = self.get_ensembl_genomes_candidate_species_folders(division_url, assembly_names)

        if not candidate_folder_urls:
            logger.error(f"[Ensembl] :: No candidate annotation folder found for {self.species}")
            return []

        annotation_files = []
        for folder_url in candidate_folder_urls:
            annotation_filename = self.get_annotation_file(folder_url)
            annotation_full_url = folder_url + annotation_filename
            logger.info(f"[Ensembl] :: Found annotation URL: {annotation_full_url}.\nDownloading...")
            annotation_file = self.target_folder / annotation_filename
            self.download_file(annotation_full_url, annotation_file)
            annotation_files.append(annotation_file)

        return annotation_files


    def get_ensembl_organisms_annotations(self) -> list[Path]:
        annotation_files = []
        species_folder_urls = self.get_ensembl_organisms_species_folders()
        for url in species_folder_urls:
            candidate_folder_urls = self.get_ensembl_organisms_subfolder(url)
            for candidate_folder_url in candidate_folder_urls:
                annotation_url = self.get_ensembl_organism_annotation_url(candidate_folder_url)
                if annotation_url is not None:
                    assembly_accession = candidate_folder_url.rstrip('/').split('/')[-1]
                    annotation_file = self.target_folder / f"{assembly_accession}.{annotation_url.split('/')[-1]}"
                    self.download_file(annotation_url, annotation_file)
                    annotation_files.append(annotation_file)
        return annotation_files


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


    def get_species_division(self, species_info: list[dict]) -> str:
        found_divisions = list({d["division"] for d in species_info})
        # this should not happen (and if it does, it's an issue on Ensembl's side)
        if not found_divisions:
            raise ValueError(f"Could not find any division for species {self.species_taxid}...")
        # we should never have multiple possible divisions for a single species
        # it is like if a species belonged to multiple kingdoms at the same time...
        if len(found_divisions) > 1:
            logger.warning(f"Multiple divisions found for species Taxon ID {self.species_taxid}: {found_divisions}. Taking the first one.")
        # there should be only one division
        return found_divisions[0]


    def get_assembly_names(self, species_info: list[dict]) -> list[str]:
        # taking all assembly names
        assembly_names = list({d.get("name", "") for d in species_info})
        return assembly_names


    def get_division_url(self, division: str) -> str:
        """
        In Ensembl, species are separated into divisions.
        Returns the URL for the division of the given species.
        """
        # the URL is different for vertebrates
        if division == self.ENSEMBL_VERTEBRATE_DIVISION:
            return self.ENSEMBL_VERTEBRATES_BASE_URL
        else:
            division_folder = self.ENSEMBL_DIVISION_TO_FOLDER[division]
            return self.ENSEMBL_GENOMES_BASE_URL.format(division_folder)


    def get_ensembl_genomes_candidate_species_folders(self, url: str, assembly_names: list[str], first_level: bool = True) -> list[str]:
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
                if folder_name.startswith(self.species) or folder_name in assembly_names:
                    folder_urls.append(folder_url)
                elif folder_name.endswith("_collection/"):
                    collection_folder_urls += self.get_ensembl_genomes_candidate_species_folders(folder_url, assembly_names, first_level=False)
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


    def get_ensembl_organisms_formatted_species(self) -> str:
        splitted = self.species.replace(" ", "_").split("_")
        return f"{splitted[0].capitalize()}_{splitted[1].lower()}"


    def get_ensembl_organisms_species_folders(self):
        formated_species = self.get_ensembl_organisms_formatted_species()
        soup = parse_page_data(self.ENSEMBL_ORGANISMS_BASE_URL)
        folder_names = []
        for item in soup.find_all("tr"):
            # all line sections
            line_sections = list(item.find_all("td"))
            if len(line_sections) < 5:
                continue
            folder_name = line_sections[1].text.strip()
            if folder_name in self.ENSEMBL_ORGANISMS_EXCLUDED_FOLDERS:
                continue
            if folder_name.startswith(formated_species):
                folder_names.append(folder_name)
        return [f"{self.ENSEMBL_ORGANISMS_BASE_URL}{folder_name}" for folder_name in folder_names]


    def get_ensembl_organisms_subfolder(self, url: str, accepted_folders: list[str] = [], excluded_folders: list[str] = []) -> list[str]:
        soup = parse_page_data(url)
        folder_names = []
        for item in soup.find_all("tr"):
            # all line sections
            line_sections = list(item.find_all("td"))
            if len(line_sections) < 5:
                continue
            folder_name = line_sections[1].text.strip()
            if accepted_folders and folder_name not in accepted_folders:
                continue
            if folder_name in excluded_folders + self.ENSEMBL_ORGANISMS_EXCLUDED_FOLDERS:
                continue
            folder_names.append(folder_name)
        return [url + folder_name for folder_name in folder_names]


    def get_ensembl_organism_annotation_url(self, url: str) -> str | None:
        first_level_folder_urls = self.get_ensembl_organisms_subfolder(url, excluded_folders=["genome", "vep"])
        if not first_level_folder_urls:
            return None
        second_level_folder_url = None
        for first_level_folder_url in first_level_folder_urls:
            second_level_folder_urls = self.get_ensembl_organisms_subfolder(first_level_folder_url, accepted_folders=["geneset/"])
            if second_level_folder_urls:
                second_level_folder_url = second_level_folder_urls[0]
                break
        if second_level_folder_url is None:
            return None
        third_level_folder_urls = self.get_ensembl_organisms_subfolder(second_level_folder_url)
        if not third_level_folder_urls:
            return None
        return third_level_folder_urls[0] + self.ENSEMBL_ORGANISMS_ANNOTATION_FILENAME


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
