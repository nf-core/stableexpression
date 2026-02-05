#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import shutil
import sys
import zipfile
from pathlib import Path

import requests
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

# Modern NCBI API
NCBI_DATASET_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/"

NCBI_TAXONOMY_ENDPOINT = "taxonomy"
NCBI_GENOME_DATASET_REPORT_BASE_ENDPOINT = "genome/taxon/{taxid}/dataset_report"
NCBI_DOWNLOAD_ENDPOINT = "genome/download"


NCBI_GENOME_DATASET_REPORT_API_PARAMS = {
    "filters.has_annotation": True,
    "page_size": 1000,
}
NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}

DOWNLOADED_FILENAME = "ncbi_dataset.zip"
ACCESSION_FILE = "accession.txt"


#####################################################
#####################################################
# PARSER
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get best assembly for a specific taxon ID"
    )
    parser.add_argument("--species", type=str, required=True, help="Species name")
    return parser.parse_args()


#####################################################
#####################################################
# REQUESTS
#####################################################
#####################################################


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_post_request_to_ncbi_dataset(endpoint: str, data: dict, params: dict = {}):
    url = NCBI_DATASET_API_URL + endpoint
    response = requests.post(url, headers=NCBI_API_HEADERS, json=data, params=params)
    response.raise_for_status()
    return response.json()


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_get_request_to_ncbi_dataset(endpoint: str, params: dict = {}):
    url = NCBI_DATASET_API_URL + endpoint
    response = requests.get(url, headers=NCBI_API_HEADERS, params=params)
    response.raise_for_status()
    return response.json()


#####################################################
#####################################################
# DATA HANDLING
#####################################################
#####################################################


def get_species_taxid(species: str) -> int:
    data = {"taxons": [species]}
    result = send_post_request_to_ncbi_dataset(NCBI_TAXONOMY_ENDPOINT, data)

    if len(result["taxonomy_nodes"]) > 1:
        raise ValueError(f"Multiple taxids for species {species}")
    metadata = result["taxonomy_nodes"][0]

    if "taxonomy" not in metadata:
        logger.info(f"Could not find taxonomy results for species {species}")
        if "errors" in metadata:
            for error in metadata["errors"]:
                logger.error(f"Error: {error['reason']}\n")
                sys.exit(100)
    return int(metadata["taxonomy"]["tax_id"])


def get_assembly_reports(taxid: int):
    result = send_get_request_to_ncbi_dataset(
        endpoint=NCBI_GENOME_DATASET_REPORT_BASE_ENDPOINT.format(taxid=taxid),
        params=NCBI_GENOME_DATASET_REPORT_API_PARAMS,
    )
    return result.get("reports", [])


def get_assembly_with_best_stats(reports: list[dict]):
    sorted_reports = sorted(
        reports,
        key=lambda x: (
            int(x.get("assembly_stats").get("total_sequence_length", 0)),
            -int(x.get("assembly_stats", {}).get("total_number_of_chromosomes", 1e9)),
        ),
        reverse=True,
    )
    return sorted_reports[0]


def get_current_assemblies(reports: list[dict]) -> dict | None:
    current_assembly_reports = [
        report
        for report in reports
        if report.get("assembly_info", {}).get("refseq_category") == "reference genome"
    ]
    if not current_assembly_reports:
        return None

    refseq_reports = [
        report
        for report in current_assembly_reports
        if report.get("source_database") == "SOURCE_DATABASE_REFSEQ"
    ]

    if refseq_reports:
        return refseq_reports[0]
    else:
        return None


def get_reference_assembly(reports: list[dict]) -> dict:
    best_assembly_report = get_current_assemblies(reports)
    if best_assembly_report is not None:
        return best_assembly_report
    else:
        return get_assembly_with_best_stats(reports)


def format_species_name(species: str):
    return species.replace("_", " ").lower()


def download_genome_annotation(genome_accession: str) -> str:
    data = {"accessions": [genome_accession], "include_annotation_type": ["GENOME_GFF"]}
    params = {"filename": DOWNLOADED_FILENAME}
    send_post_request_to_ncbi_dataset(NCBI_TAXONOMY_ENDPOINT, data, params)


def extract_annotation_file_from_archive():
    with zipfile.ZipFile(DOWNLOADED_FILENAME, "r") as zip_ref:
        zip_ref.extractall()

    valid_files = list(Path().cwd().glob(f"ncbi_dataset/data/{accession}/*.gff"))

    if not valid_files:
        raise ValueError(f"No annotation file found for accession {accession}")

    if len(valid_files) > 1:
        logger.warning(
            f"Multiple annotation files found for accession {accession}. Taking the first one"
        )

    annotation_file = valid_files[0]
    shutil.move(annotation_file, f"{accession}.gff")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    species = format_species_name(args.species)

    species_taxid = get_species_taxid(species)
    logger.info(f"Species taxid: {species_taxid}")

    logger.info(f"Getting best NCBI assembly for taxid: {species_taxid}")
    reports = get_assembly_reports(species_taxid)

    if not reports:
        logger.error(f"No assembly reports found for taxid {species_taxid}")
        sys.exit(100)

    # looping while we can get an annotation file
    annotation_found = False
    while not annotation_found:
        best_assembly_report = get_reference_assembly(reports)
        logger.info(
            f"Best assembly: {best_assembly_report['accession']}. Trying to download annotation"
        )
        accession = best_assembly_report["accession"]
        try:
            download_genome_annotation(accession)
            extract_annotation_file_from_archive()
            annotation_found = True
        except Exception as e:
            logger.error(f"Error downloading annotation for accession {accession}: {e}")

        if not annotation_found:
            # Remove the best assembly report from the list of reports
            reports = [report for report in reports if report["accession"] != accession]

    if not annotation_found:
        logger.error(f"No annotation found for taxid {species_taxid}")
        sys.exit(100)

    logger.info("Done")
