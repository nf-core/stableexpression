#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
import shutil
import zipfile
from pathlib import Path

import httpx
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logger = logging.getLogger(__name__)
logging.getLogger("httpx").setLevel(logging.ERROR)

NCBI_DATASET_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
NCBI_BASE_HEADERS = {"accept": "application/json", "content-type": "application/json"}

NCBI_TAXONOMY_ENDPOINT = "taxonomy"

NCBI_GENOME_REPORT_BASE_ENDPOINT = "genome/taxon/{taxid}/dataset_report"
NCBI_GENOME_REPORT_API_PARAMS = {
    "filters.has_annotation": True,
    "page_size": 1000,
}

NCBI_GENOME_DOWNLOAD_ENDPOINT = "genome/download"
NCBI_GENOME_DOWNLOAD_HEADERS = {
    "accept": "application/zip",
    "content-type": "application/json",
}
NCBI_GENOME_DOWNLOAD_BASE_DATA = {
    "include_annotation_type": ["GENOME_GFF", "GENOME_GTF"]
}
DOWNLOADED_FILENAME = "ncbi_dataset.zip"

STOP_RETRY_AFTER_DELAY = 120


#####################################################
#####################################################
# httpx
#####################################################
#####################################################


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(data: dict):
    url = NCBI_DATASET_API_URL + NCBI_TAXONOMY_ENDPOINT
    with httpx.Client() as client:
        response = client.post(url, headers=NCBI_BASE_HEADERS, json=data)
        response.raise_for_status()
        return response.json()


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_genome_report(taxid: str | int):
    url = NCBI_DATASET_API_URL + NCBI_GENOME_REPORT_BASE_ENDPOINT.format(
        taxid=str(taxid)
    )
    with httpx.Client() as client:
        response = client.get(url, headers=NCBI_BASE_HEADERS)
        response.raise_for_status()
        return response.json()


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_genome_download(genome_accession: str):
    url = NCBI_DATASET_API_URL + NCBI_GENOME_DOWNLOAD_ENDPOINT
    data = NCBI_GENOME_DOWNLOAD_BASE_DATA | {"accessions": [genome_accession]}
    with httpx.Client() as client:
        response = client.post(
            url,
            headers=NCBI_BASE_HEADERS,
            json=data
        )
        response.raise_for_status()
        # write response content to file
        outfile = Path(DOWNLOADED_FILENAME)
        with open(outfile, "wb") as f:
            f.write(response.content)
    return outfile


#####################################################
#####################################################
# DATA HANDLING
#####################################################
#####################################################


class NCBIAnnotationManager:

    species: str
    species_taxid: int
    reports: list[dict]

    def __init__(self, species: str):
        self.species = species
        self.species_taxid = self.get_species_taxid(self.species)
        logger.info(f"[NCBI] :: Getting reference genome reports for taxid: {self.species_taxid}")
        self.reports = self.get_assembly_reports()


    @staticmethod
    def get_species_taxid(species: str) -> int:
        formatted_species = species.replace("_", " ").lower()
        data = {"taxons": [formatted_species]}
        result = send_request_to_ncbi_taxonomy(data)
        if len(result["taxonomy_nodes"]) > 1:
            raise ValueError(f"Multiple taxids for species {species}")
        metadata = result["taxonomy_nodes"][0]
        if "taxonomy" not in metadata:
            raise ValueError(f"Could not find taxonomy results for species {species}")
        return int(metadata["taxonomy"]["tax_id"])


    def get_assembly_reports(self):
        result = send_request_to_ncbi_genome_report(self.species_taxid)
        return result.get("reports", [])


    def get_sorted_reference_genome_reports(self, refseq_only: bool = False) -> list[dict]:
        # selecting genome annotated as 'reference'
        reference_reports = [
            report
            for report in self.reports
            if report.get("assembly_info", {}).get("refseq_category") == "reference genome"
        ]
        if refseq_only:
            reference_reports = [
                report
                for report in reference_reports
                if report.get("source_database") == "SOURCE_DATABASE_REFSEQ"
            ]
        # sorting by assembly statistics
        return sorted(
            reference_reports,
            key=lambda x: (
                int(x.get("assembly_stats", {}).get("total_sequence_length", 0)),
                -int(x.get("assembly_stats", {}).get("total_number_of_chromosomes", 1e9)),
            ),
            reverse=True,
        )


    @staticmethod
    def download_genome_annotation(genome_accession: str) -> Path:
        download_archive = send_request_to_ncbi_genome_download(genome_accession)
        if not download_archive.exists():
            raise FileNotFoundError(
                f"Downloaded file not found for accession {genome_accession}"
            )
        return download_archive


    @staticmethod
    def extract_annotation_file_from_archive(archive: Path, accession: str, folder: Path) -> Path:
        logger.info(f"Extracting annotation file from archive {archive}")
        with zipfile.ZipFile(archive, "r") as zip_ref:
            for file in zip_ref.namelist():
                if file.endswith(('.gff', '.gtf')):
                    zip_ref.extract(file, path='.')

        valid_files = list(Path().cwd().glob(f"ncbi_dataset/data/{accession}/*.gff"))
        if not valid_files:
            valid_files = list(Path().cwd().glob(f"ncbi_dataset/data/{accession}/*.gtf"))
            if not valid_files:
                raise ValueError(f"No annotation file found for accession {accession}")

        if len(valid_files) > 1:
            logger.warning(
                f"Multiple annotation files found for accession {accession}. Taking the first one"
            )

        folder.mkdir(parents=True, exist_ok=True)
        annotation_file = folder / f"{accession}.gff"
        downloaded_file = valid_files[0]
        shutil.move(downloaded_file, annotation_file)

        logger.info(f"Annotation file for accession {accession} saved to {annotation_file}. Removing archive and extracted folder")
        shutil.rmtree("ncbi_dataset")
        archive.unlink()

        return annotation_file
