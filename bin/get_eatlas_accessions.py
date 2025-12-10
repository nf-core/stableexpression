#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import random
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

ALLOWED_PLATFORMS = ["rnaseq", "microarray"]
# accessions that should not be fetched automatically:
# - E-GTEX-8 contains 17350 samples (way too big)
EXCLUDED_ACCESSION_PATTERNS = ["E-GTEX-"]

ALL_EXP_URL = "https://www.ebi.ac.uk/gxa/json/experiments/"
ACCESSION_OUTFILE_NAME = "accessions.txt"
# ALL_EXPERIMENTS_METADATA_OUTFILE_NAME = "all_experiments.metadata.tsv"
SPECIES_EXPERIMENTS_METADATA_OUTFILE_NAME = "species_experiments.metadata.tsv"
SELECTED_EXPERIMENTS_METADATA_OUTFILE_NAME = "selected_experiments.metadata.tsv"
FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME = "filtered_experiments.keywords.yaml"

SAMPLING_QUOTA_OUTFILE = "sampling_quota.txt"


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
    parser.add_argument(
        "--platform", type=str, help="Platform type", choices=ALLOWED_PLATFORMS
    )
    parser.add_argument(
        "--random-sampling-size",
        dest="random_sampling_size",
        type=int,
        help="Random sampling size",
    )
    parser.add_argument(
        "--random-sampling-seed",
        dest="random_sampling_seed",
        type=int,
        help="Random sampling seed",
    )
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


def filter_by_platform(experiments: list[dict], platform: str | None):
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
            # parsed_platform is in ["rnaseq", "microarray", "proteomics", ...]
            parsed_platform = (
                parsed_technology_type.lower().split(" ")[0].replace("-", "")
            )

            if platform is not None:
                if parsed_platform == platform:
                    platform_experiments.append(exp_dict)
            else:
                if parsed_platform in ALLOWED_PLATFORMS:
                    platform_experiments.append(exp_dict)

        else:
            logger.warning(
                f"Technology type not found for experiment {exp_dict['accession']}"
            )
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


def filter_out_excluded_accessions(experiments: list[dict]) -> list[dict]:
    valid_experiments = []
    for exp_dict in experiments:
        for accession_pattern in EXCLUDED_ACCESSION_PATTERNS:
            if exp_dict["experimentAccession"].startswith(accession_pattern):
                logger.warning(
                    f"Skipping experiment {exp_dict['experimentAccession']} due to exclusion pattern"
                )
                break
        else:
            valid_experiments.append(exp_dict)
    return valid_experiments


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


def sample_experiments_randomly(
    experiments: list[dict], sampling_size: int, seed: int
) -> tuple[list[str], bool]:
    random.seed(seed)
    sampled_experiments = []

    total_nb_samples = 0
    sampling_quota_reached = False
    experiments_left = list(experiments)
    while experiments_left:
        # if the min number of samples is greater than the remaining space left, we get out of the loop
        experiments_left_nb_samples = [exp["nb_samples"] for exp in experiments_left]
        min_nb_samples = min(experiments_left_nb_samples)
        if min_nb_samples > sampling_size - total_nb_samples:
            sampling_quota_reached = True
            logger.warning("Sampling quota reached")
            break

        experiment = None
        test_total_nb_samples = int(total_nb_samples)
        experiments_not_tested = list(experiments_left)
        while experiments_not_tested:
            experiment = random.choice(experiments_not_tested)
            experiments_not_tested.remove(experiment)
            # if we do not exceed the sampling size with this experiment
            # we keep it
            test_total_nb_samples = total_nb_samples + experiment["nb_samples"]
            if test_total_nb_samples <= sampling_size:
                break

        # this should not happen but we keep it for safety
        if experiment is None:
            logger.error("No experiment found")
            continue

        total_nb_samples = test_total_nb_samples
        experiments_left.remove(experiment)
        sampled_experiments.append(experiment)

    return [exp["accession"] for exp in sampled_experiments], sampling_quota_reached


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
    experiments = get_eatlas_experiments()

    logger.info("Filtering on species name")
    experiments = get_species_experiments(experiments, species_name)
    logger.info(f"Found {len(experiments)} experiments for species {species_name}")

    logger.info("Filtering experiments based on platform")
    experiments = filter_by_platform(experiments, args.platform)

    logger.info("Filtering out excluded accessions")
    experiments = filter_out_excluded_accessions(experiments)

    logger.info("Parsing experiments")
    with Pool(processes=args.nb_cpus) as pool:
        results = pool.map(parse_experiment, experiments)

    if keywords:
        logger.info(f"Filtering experiments with keywords {keywords}")
        func = partial(filter_experiment_with_keywords, keywords=keywords)
        with Pool(processes=args.nb_cpus) as pool:
            results = [res for res in pool.map(func, results) if res is not None]
        logger.info(
            f"Found {len(results)} experiments corresponding to keywords {keywords}"
        )

    # getting accessions of selected experiments
    selected_accessions = [exp_dict["accession"] for exp_dict in results]

    if args.random_sampling_size and args.random_sampling_seed:
        selected_accession_to_nb_samples = [
            {
                "accession": exp_dict["experimentAccession"],
                "nb_samples": exp_dict["numberOfAssays"],
            }
            for exp_dict in experiments
            if exp_dict["experimentAccession"] in selected_accessions
        ]

        nb_samples_df = pd.DataFrame.from_dict(selected_accession_to_nb_samples)
        nb_samples_df.to_csv("selected_accession_to_nb_samples.csv", index=False)

        logger.info("Sampling experiments randomly")
        selected_accessions, sampling_quota_reached = sample_experiments_randomly(
            selected_accession_to_nb_samples,
            args.random_sampling_size,
            args.random_sampling_seed,
        )
        logger.info(
            f"Kept {len(selected_accessions)} experiments after random sampling"
        )

        # writing status to file
        # so that the wrapper module can get the status
        with open(SAMPLING_QUOTA_OUTFILE, "w") as fout:
            sampling_status = "full" if sampling_quota_reached else "ok"
            fout.write(sampling_status)

    # keeping metadata only for selected experiments
    selected_experiments = get_metadata_for_selected_experiments(experiments, results)

    if not selected_accessions:
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
    df = pd.DataFrame.from_dict(experiments)
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

    if results:
        # exporting list of selected experiments with their keywords
        logger.info(
            f"Writing filtered experiments with keywords to {FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME}"
        )
        with open(FILTERED_EXPERIMENTS_WITH_KEYWORDS_OUTFILE_NAME, "w") as fout:
            yaml.dump(results, fout)


if __name__ == "__main__":
    main()
