#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys

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
NCBI_TAXONOMY_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
NCBI_GENOME_DATASET_REPORT_API_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report"
)
NCBI_GENOME_DATASET_REPORT_API_PARAMS = "filters.has_annotation=true&page_size=1000"
NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}

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
def send_request_to_ncbi_taxonomy(taxid: str | int):
    taxons = [str(taxid)]
    data = {"taxons": taxons}
    response = requests.post(NCBI_TAXONOMY_API_URL, headers=NCBI_API_HEADERS, json=data)
    response.raise_for_status()
    return response.json()


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_genome_dataset_report_api(taxid: int):
    url = NCBI_GENOME_DATASET_REPORT_API_URL.format(taxid=taxid)
    url += f"?{NCBI_GENOME_DATASET_REPORT_API_PARAMS}"
    response = requests.get(url, headers=NCBI_API_HEADERS)
    response.raise_for_status()
    return response.json()


#####################################################
#####################################################
# DATA HANDLING
#####################################################
#####################################################


def get_species_taxid(species: str) -> int:
    result = send_request_to_ncbi_taxonomy(species)

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
    result = send_request_to_ncbi_genome_dataset_report_api(species_taxid)

    try:
        reports = result["reports"]
        best_assembly_report = get_reference_assembly(reports)
        logger.info(f"Best assembly: {best_assembly_report['accession']}")
    except Exception as e:
        logger.error(f"Could not get any assembly for taxid {species_taxid}: {e}")
        sys.exit(100)

    with open(ACCESSION_FILE, "w") as fout:
        fout.write(best_assembly_report["accession"])

    logger.info("Done")
