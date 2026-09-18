#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
import shutil
import pandas as pd
from tenacity import RetryError
from ensembl_annotation_manager import EnsemblAnnotationManager
from ncbi_annotation_manager import NCBIAnnotationManager

from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL_ANNOTATION_LOCAL_FOLDER = "ensembl_annotations"
NCBI_ANNOTATION_LOCAL_FOLDER = "ncbi_annotations"

GFF3_COLUMN_DTYPES = {
    "chromosome": str,
    "source": str,
    "feature": str,
    "start": int,
    "end": int,
    "score": str,
    "strand": str,
    "phase": str,
    "attributes": str,
}

OUTPUT_DIR = "selected"


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
    parser.add_argument(
        "--gene-ids",
        type=Path,
        dest="gene_ids_file",
        required=True,
        help="File containing gene IDs",
    )
    parser.add_argument(
        "--skip-ensembl",
        action="store_true",
        dest="skip_ensembl"
    )
    parser.add_argument(
        "--skip-ncbi",
        action="store_true",
        dest="skip_ncbi"
    )
    return parser.parse_args()


def get_ensembl_annotations(species: str) -> list[Path]:

    logger.info("[Ensembl] :: Searching for the appropriate annotation folder")
    candidate_folder_urls = EnsemblAnnotationManager(species).get_candidate_folders()

    if not candidate_folder_urls:
        logger.error(f"[Ensembl] :: No candidate annotation folder found for {species}")
        return []

    target_folder = Path(ENSEMBL_ANNOTATION_LOCAL_FOLDER)
    target_folder.mkdir(parents=True, exist_ok=True)

    annotation_files = []
    for folder_url in candidate_folder_urls:
        annotation_filename = EnsemblAnnotationManager.get_annotation_file(folder_url)

        annotation_full_url = folder_url + annotation_filename
        logger.info(f"[Ensembl] :: Found annotation URL: {annotation_full_url}.\nDownloading...")
        annotation_file = target_folder / annotation_filename
        EnsemblAnnotationManager.download_file(annotation_full_url, annotation_file)
        annotation_files.append(annotation_file)

    return annotation_files


def get_ncbi_annotations(species: str) -> list[Path]:
    """
    Downloading genome annotation
    """

    target_folder = Path(NCBI_ANNOTATION_LOCAL_FOLDER)
    target_folder.mkdir(parents=True, exist_ok=True)

    # getting list of assemblies from the most 'reference ' to the least
    reference_reports = NCBIAnnotationManager(species).get_sorted_reference_genome_reports()

    if not reference_reports:
        logger.info(f"[NCBI] :: No assembly reports found for species {species}")
        return []

    annotation_files = []
    for report in reference_reports:
        accession = report['accession']
        logger.info(f"[NCBI] :: Trying to download annotation of assembly {accession}.")
        try:
            download_archive = NCBIAnnotationManager.download_genome_annotation(accession)
            annotation_file = NCBIAnnotationManager.extract_annotation_file_from_archive(download_archive, accession, target_folder)
            annotation_files.append(annotation_file)
        except Exception as e:
            logger.error(f"[NCBI] :: Error downloading annotation for accession {accession}: {e}")

    if not annotation_files:
        logger.error(f"[NCBI] :: No annotation found for species {species}")

    return annotation_files


def parse_gene_ids(file: Path) -> list[str]:
    with open(file, "r") as fin:
        return sorted({line.strip() for line in fin})


def parse_annotation(annotation_file: Path) -> pd.DataFrame:
    return pd.read_csv(
        annotation_file,
        sep="\t",
        names=list(GFF3_COLUMN_DTYPES.keys()),
        dtype=GFF3_COLUMN_DTYPES,
        comment="#",
        on_bad_lines="warn",
    )


def parse_gene_ids_from_annotation(file: Path) -> list[str]:
    """
    Extract gene ID from attributes column for each gene feature
    """
    logger.info(f"Parsing gene IDs from annotation file {file}")
    df = parse_annotation(file)

    suffix = file.suffix if file.suffix != '.gz' else file.suffixes[-2]
    if suffix in [".gff3", ".gff"]:
        pattern = r"ID=gene[:\-]([^;]+)"
    elif suffix == ".gtf":
        pattern = r'gene_id\s+"([^"]*)"'
    else:
        raise ValueError(f"Unsupported file suffix: {suffix}")

    return (
        df.loc[df["feature"] == 'gene']['attributes']
        .str.extract(pattern, expand=False)
        .drop_duplicates()
        .dropna()
        .tolist()
    )


def get_annotation_matching_gene_ids(annotations: list[Path], unique_gene_ids: list[str]) -> Path | None:
    """
    Find the annotation file that contains the most gene IDs in common with the list provided.
    """
    annotation_file_to_nb_common_gene_ids = {}
    for annotation_file in annotations:

        annotation_gene_ids = parse_gene_ids_from_annotation(annotation_file)
        if annotation_gene_ids:
            genes_to_show = [str(gene_id) for gene_id in annotation_gene_ids[:min(3, len(annotation_gene_ids))]]
            logger.info(f"{annotation_file.name} :: found {len(annotation_gene_ids)} gene IDs like {', '.join(genes_to_show)}")

            gene_id_intersection = set(annotation_gene_ids).intersection(unique_gene_ids)
            logger.info(f"{annotation_file.name} :: {len(gene_id_intersection)} gene IDs in commmon")

            annotation_file_to_nb_common_gene_ids[annotation_file] = len(gene_id_intersection)

    if not annotation_file_to_nb_common_gene_ids:
        logger.error("Could not parse any gene IDs from downloaded annotations...")
        return None

    max_nb_common_gene_ids = max(annotation_file_to_nb_common_gene_ids.values())
    if max_nb_common_gene_ids == 0:
        logger.error("Could not find any annotation file having gene IDs in common with list provided...")
        return None

    # if multiple annotation file share the same number of common gene IDs
    # taking the first one
    best_annotations = [
        file for file, nb_common_gene_ids in annotation_file_to_nb_common_gene_ids.items()
        if nb_common_gene_ids == max_nb_common_gene_ids
    ]
    return best_annotations[0]


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    args = parse_args()

    species = args.species
    unique_gene_ids = parse_gene_ids(args.gene_ids_file)
    logger.info(f"Got {len(unique_gene_ids)} unique gene IDs like {', '.join(unique_gene_ids[:3])}")

    try:
        ##################################################################
        # ENSEMBL
        ##################################################################

        selected_annotation = None
        search_on_ncbi = True

        if not args.skip_ensembl:
            logger.info(f"[Ensembl] :: Fetching annotations for {species}")
            ensembl_annotations = get_ensembl_annotations(species)
            selected_annotation = get_annotation_matching_gene_ids(
                ensembl_annotations,
                unique_gene_ids
            )
            if selected_annotation is not None:
                logger.warning("Could not find any suitable annotation in Ensembl. Trying with NCBI")
                search_on_ncbi = False

        ##################################################################
        # NCBI
        ##################################################################

        if search_on_ncbi and not args.skip_ncbi:
            logger.info(f"[NCBI] :: Fetching annotations for {species}")
            ncbi_annotations = get_ncbi_annotations(species)
            selected_annotation = get_annotation_matching_gene_ids(
                ncbi_annotations,
                unique_gene_ids
            )

    except RetryError as e: # if servers are unavailable for some reason
        logger.exception(e)
        sys.exit(101) # retrying at the module level


    ## ################################################################
    # MOVING THE SELECTED ANNOTATION TO THE PLACE
    ##################################################################

    if selected_annotation is None:
        raise ValueError(f"Could not find any annotation for species {species}.")

    logger.info(f"Selected annotation file: {selected_annotation.name}")
    outdir = Path().cwd() / OUTPUT_DIR
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / selected_annotation.name
    shutil.move(selected_annotation, outfile)

    logger.info("Done")


if __name__ == "__main__":
    main()
