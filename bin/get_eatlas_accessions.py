#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from functools import partial
from multiprocessing import Pool

import pandas as pd
import requests
import yaml
from natural_language_utils import keywords_in_fields
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ALL_EXP_URL = "https://www.ebi.ac.uk/gxa/json/experiments/"
ACCESSION_OUTFILE_NAME = "accessions.txt"
# ALL_EXPERIMENTS_METADATA_OUTFILE_NAME = "all_experiments.metadata.tsv"
SPECIES_EXPERIMENTS_METADATA_OUTFILE_NAME = "species_experiments.metadata.tsv"
SELECTED_EXPERIMENTS_METADATA_OUTFILE_NAME = "selected_experiments.metadata.tsv"
FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME = "filtered_experiments.keywords.yaml"


##################################################################
##################################################################
# FUNCTIONS
##################################################################
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Get expression atlas accessions")
    parser.add_argument(
        "--species",
        type=str,
        required=True,
        help="Search Expression Atlas for this specific species",
    )
    parser.add_argument(
        "--cpus",
        dest="nb_cpus",
        type=int,
        required=True,
        help="Number of CPUs to use",
    )
    parser.add_argument(
        "--keywords",
        type=str,
        nargs="*",
        help="Keywords to search for in experiment description",
    )
    parser.add_argument("--platform", type=str, help="Platform type")
    return parser.parse_args()


@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def get_data(url: str) -> dict:
    """
    Queries a URL and returns the data as a JSON object

    Parameters
    ----------
    url : str
        The URL to query

    Returns
    -------
    data : dict
        The JSON object returned by the query

    Raises
    ------
    RuntimeError
        If the query fails
    """
    response = requests.get(url)
    response.raise_for_status()
    return response.json()


def get_experiment_description(exp_dict: dict):
    """
    Gets the description from an experiment dictionary

    Parameters
    ----------
    exp_dict : dict
        The experiment dictionary

    Returns
    -------
    description : str
        The experiment description

    Raises
    ------
    KeyError
        If the description field is not found in the experiment dictionary
    """
    if "experiment" in exp_dict:
        if "description" in exp_dict["experiment"]:
            return exp_dict["experiment"]["description"]
        else:
            raise KeyError(f"Could not find description field in {exp_dict}")
    elif "experimentDescription" in exp_dict:
        return exp_dict["experimentDescription"]
    else:
        raise KeyError(f"Could not find description field in {exp_dict}")


def get_experiment_accession(exp_dict: dict):
    """
    Gets the accession from an experiment dictionary

    Parameters
    ----------
    exp_dict : dict
        The experiment dictionary

    Returns
    -------
    accession : str
        The experiment accession

    Raises
    ------
    KeyError
        If the accession field is not found in the experiment dictionary
    """
    if "experiment" in exp_dict:
        if "accession" in exp_dict["experiment"]:
            return exp_dict["experiment"]["accession"]
        else:
            raise KeyError(f"Could not find accession field in {exp_dict}")
    elif "experimentAccession" in exp_dict:
        return exp_dict["experimentAccession"]
    else:
        raise KeyError(f"Could not find accession field in {exp_dict}")


def get_properties_values(exp_dict: dict):
    """
    Gets all values from properties from an experiment dictionary

    Parameters
    ----------
    exp_dict : dict
        The experiment dictionary

    Returns
    -------
    values : list
        A list of all values from properties
    """
    values = []
    for column_header_dict in exp_dict["columnHeaders"]:
        key_found = False
        for key in ["assayGroupSummary", "contrastSummary"]:
            if key in column_header_dict:
                for property_dict in column_header_dict[key]["properties"]:
                    values.append(property_dict["testValue"])
                key_found = True
                break
        if not key_found:
            raise KeyError(f"Could not find property value in {column_header_dict}")
    # removing empty strings
    values = [value for value in values if value != ""]
    # removing duplicates
    return list(set(values))


def get_eatlas_experiments():
    """
    Gets all experiments from Expression Atlas

    Parameters
    ----------

    Returns
    -------
    experiments : list
        A list of experiment dictionaries
    """
    data = get_data(ALL_EXP_URL)
    return data["experiments"]


def get_platform_specific_experiments(experiments: list[dict], platform: str):
    """
    Gets all experiments for a given platform from Expression Atlas
    Possible platforms in Expression Atlas are 'rnaseq', 'microarray', 'proteomics'

    Parameters
    ----------
    experiments: list[str]
    platform : str
        Name of platform. Example: "rnaseq"

    Returns
    -------
    experiments : list
        A list of experiment dictionaries
    """
    platform_experiments = []
    for exp_dict in experiments:
        if technology_type := exp_dict.get("technologyType"):
            parsed_technology_type = (
                technology_type[0]
                if isinstance(technology_type, list)
                else technology_type
            )
            parsed_platform = (
                parsed_technology_type.lower().split(" ")[0].replace("-", "")
            )
            if platform == parsed_platform:
                platform_experiments.append(exp_dict)
    return platform_experiments


def get_species_experiments(experiments: list[dict], species: str):
    """
    Gets all experiments for a given species from Expression Atlas

    Parameters
    ----------
    experiments: list[str]
    species : str
        Name of species. Example: "Arabidopsis thaliana"

    Returns
    -------
    experiments : list
        A list of experiment dictionaries
    """
    species_experiments = []
    for exp_dict in experiments:
        if exp_dict["species"] == species:
            species_experiments.append(exp_dict)
    return species_experiments


def get_experiment_data(exp_dict: dict):
    """
    Gets the full data for an experiment given its dictionary

    Parameters
    ----------
    exp_dict : dict
        The experiment dictionary

    Returns
    -------
    exp_data : dict
        The full experiment data
    """
    exp_url = ALL_EXP_URL + exp_dict["experimentAccession"]
    return get_data(exp_url)


def parse_experiment(exp_dict: dict):
    # getting accession and description
    accession = get_experiment_accession(exp_dict)
    description = get_experiment_description(exp_dict)
    # getting properties of this experiment
    exp_data = get_experiment_data(exp_dict)
    properties_values = get_properties_values(exp_data)

    return {
        "accession": accession,
        "description": description,
        "properties": properties_values,
    }


def filter_experiment_with_keywords(exp_dict: dict, keywords: list[str]) -> dict | None:
    all_searchable_fields = [exp_dict["description"]] + exp_dict["properties"]
    found_keywords = keywords_in_fields(all_searchable_fields, keywords)
    # only returning experiments if found keywords
    if found_keywords:
        exp_dict["found_keywords"] = list(set(found_keywords))
        return exp_dict
    else:
        return None


def get_metadata_for_selected_experiments(
    experiments: list[dict], results: list[dict]
) -> list[dict]:
    filtered_accessions = [result_dict["accession"] for result_dict in results]
    return [
        exp_dict
        for exp_dict in experiments
        if get_experiment_accession(exp_dict) in filtered_accessions
    ]


def format_species_name(species: str) -> str:
    return species.replace("_", " ").capitalize().strip()


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    args = parse_args()

    results = None
    selected_accessions = []
    selected_experiments = []

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING EXPRESSION ATLAS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Getting arguments
    species_name = format_species_name(args.species)
    keywords = args.keywords

    logger.info(f"Getting experiments corresponding to species {species_name}")
    all_experiments = get_eatlas_experiments()

    if args.platform:
        logger.info(f"Getting experiments corresponding to platform {args.platform}")
        all_experiments = get_platform_specific_experiments(
            all_experiments, args.platform
        )

    species_experiments = get_species_experiments(all_experiments, species_name)
    logger.info(
        f"Found {len(species_experiments)} experiments for species {species_name}"
    )

    logger.info("Parsing experiments")
    with Pool(processes=args.nb_cpu) as pool:
        results = pool.map(parse_experiment, species_experiments)

    if keywords:
        logger.info(f"Filtering experiments with keywords {keywords}")
        func = partial(filter_experiment_with_keywords, keywords=keywords)
        with Pool(processes=args.nb_cpu) as pool:
            results = [res for res in pool.map(func, results) if res is not None]

    if results:
        logger.info(f"Kept {len(results)} experiments")
        # getting accessions of selected experiments
        selected_accessions = [exp_dict["accession"] for exp_dict in results]
        # keeping metadata only for selected experiments
        selected_experiments = get_metadata_for_selected_experiments(
            species_experiments, results
        )

    else:
        logger.warning(
            f"Could not find experiments for species {species_name} and keywords {keywords}"
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXPORTING DATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # exporting list of accessions
    logger.info(f"Writing accessions to {ACCESSION_OUTFILE_NAME}")
    with open(ACCESSION_OUTFILE_NAME, "w") as fout:
        fout.writelines([f"{acc}\n" for acc in selected_accessions])

    # exporting metadata
    logger.info(
        f"Writing metadata of all experiments for species {species_name} to {SPECIES_EXPERIMENTS_METADATA_OUTFILE_NAME}"
    )
    df = pd.DataFrame.from_dict(species_experiments)
    df.to_csv(
        SPECIES_EXPERIMENTS_METADATA_OUTFILE_NAME, sep="\t", index=False, header=True
    )

    if selected_experiments:
        logger.info(
            f"Writing metadata of filtered experiments to {SELECTED_EXPERIMENTS_METADATA_OUTFILE_NAME}"
        )
        df = pd.DataFrame.from_dict(selected_experiments)
        df.to_csv(
            SELECTED_EXPERIMENTS_METADATA_OUTFILE_NAME,
            sep="\t",
            index=False,
            header=True,
        )

    if results is not None:
        # exporting list of selected experiments with their keywords
        logger.info(
            f"Writing filtered experiments with keywords to {FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME}"
        )
        with open(FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME, "w") as fout:
            yaml.dump(results, fout)


if __name__ == "__main__":
    main()
