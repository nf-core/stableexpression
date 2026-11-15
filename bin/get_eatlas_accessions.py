#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import random
import asyncio

import httpx
import csv
from natural_language_utils import keywords_in_fields
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(level=logging.INFO)
logging.getLogger("httpx").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

# ----------------------------------------
# HTTPX
# ----------------------------------------

MAX_CONCURRENT_REQUESTS = 20
STOP_AFTER_DELAY = 60 * 30
TIMEOUT = 60 * 10 # timeout of 10 minutes since the EBI may take a while to respond

# ----------------------------------------

ALLOWED_PLATFORMS = ["rnaseq", "microarray"]
METADATA_KEYS = {
    "technologyType": "technology",
    "numberOfAssays": "mumber_of_samples",
    "experimentalFactors": "experimental_factors",
    "experimentType": "experiment_type",
    "numberOfContrasts": "number_of_contrasts",
}

ALL_EXP_URL = "https://www.ebi.ac.uk/gxa/json/experiments/"

ACCESSION_OUTFILE_NAME = "accessions.txt"
SELECTED_EXPERIMENTS_METADATA_OUTFILE_NAME = "selected_experiments.metadata.tsv"
REJECTED_EXPERIMENTS_METADATA_OUTFILE_NAME = "rejected_experiments.metadata.tsv"

SAMPLING_QUOTA_OUTFILE = "sampling_quota.txt"

NO_KEYWORD_FOUND_FAILURE_REASON = "No keyword found"
REJECTED_FOR_SUBSAMPLING_FAILURE_REASON = "Rejected due to subsampling"


# -------------------------------------------------------------------------
# FUNCTIONS
# -------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser("Get expression atlas accessions")
    parser.add_argument(
        "--species",
        type=str,
        required=True,
        help="Search Expression Atlas for this specific species",
    )
    parser.add_argument(
        "--keywords",
        type=str,
        help="Keywords to search for in experiment description (comma separated)",
        default=""
    )
    parser.add_argument(
        "--platform",
        type=str,
        help="Platform type",
        choices=ALLOWED_PLATFORMS,
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
    parser.add_argument(
        "--accessions",
        type=str,
        help="[For dev purposes / testing: provide directly accessions (separated by commas) and try to get their metadata]",
    )
    return parser.parse_args()


@retry(
    stop=stop_after_delay(STOP_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def get_eatlas_experiments() -> list[dict]:
    """Get all experiments from Expression Atlas.

    Returns
    -------
    experiments : list
        A list of experiment dictionaries

    """
    with httpx.Client(timeout=TIMEOUT) as client:
        response = client.get(ALL_EXP_URL)
        response.raise_for_status()
        return response.json()["experiments"]


@retry(
    stop=stop_after_delay(STOP_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
async def get_data(client: httpx.AsyncClient, url: str) -> dict:
    """Query a URL and returns the data as a JSON object.

    Parameters
    ----------
    client : httpx.AsyncClient
        The HTTP client to use
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
    response = await client.get(url)
    response.raise_for_status()
    return response.json()


async def parse_details_of_all_experiments(experiments: list[dict]) -> list[dict]:
    """Parse all experiments in parallel using asyncio.

    Parameters
    ----------
    experiments : list[dict]
        A list of experiment dictionaries.

    Returns
    -------
    list[dict]
        A list of parsed experiment dictionaries.

    """
    semaphore = asyncio.Semaphore(MAX_CONCURRENT_REQUESTS)
    async with httpx.AsyncClient(timeout=TIMEOUT) as client:
        tasks = [parse_experiment_details(client, exp, semaphore) for exp in experiments]
        return await asyncio.gather(*tasks)


def get_experiment_description(exp_dict: dict) -> str:
    """Get the description from an experiment dictionary.

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
        raise KeyError(f"Could not find description field in {exp_dict}")
    if "experimentDescription" in exp_dict:
        return exp_dict["experimentDescription"]
    raise KeyError(f"Could not find description field in {exp_dict}")


def get_experiment_accession(exp_dict: dict) -> str:
    """Get the accession from an experiment dictionary.

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
        raise KeyError(f"Could not find accession field in {exp_dict}")
    if "experimentAccession" in exp_dict:
        return exp_dict["experimentAccession"]
    raise KeyError(f"Could not find accession field in {exp_dict}")


def get_properties_values(exp_dict: dict) -> list:
    """Get all values from properties from an experiment dictionary.

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


def filter_by_platform(experiments: list[dict], platform: str | None) -> list[dict]:
    """Get all experiments for a given platform from Expression Atlas.

    Possible platforms in Expression Atlas are 'rnaseq', 'microarray', 'proteomics'

    Parameters
    ----------
    experiments: list[str]
        A list of experiment dictionaries
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
                technology_type[0] if isinstance(technology_type, list) else technology_type
            )
            # parsed_platform is in ["rnaseq", "microarray", "proteomics", ...]
            parsed_platform = parsed_technology_type.lower().split(" ")[0].replace("-", "")

            if platform is not None:
                if parsed_platform == platform:
                    platform_experiments.append(exp_dict)
            elif parsed_platform in ALLOWED_PLATFORMS:
                platform_experiments.append(exp_dict)

        else:
            logger.warning(
                f"Technology type not found for experiment {exp_dict['accession']}",
            )
    return platform_experiments


def get_species_name_without_subspecies(species: str) -> str:
    """Return the species name without the subspecies part.

    Ex: Hordeum vulgare subsp. vulgare -> Hordeum vulgare
    """
    return " ".join(species.split(" ")[:2])


def get_species_experiments(experiments: list[dict], species: str) -> list[dict]:
    """Get all experiments for a given species from Expression Atlas.

    Parameters
    ----------
    experiments: list[str]
        A list of experiment dictionaries
    species : str
        Name of species. Example: "Arabidopsis thaliana"

    Returns
    -------
    experiments : list
        A list of experiment dictionaries

    """
    species_experiments = []
    for exp_dict in experiments:
        # in case the species name contains a subspecies part, we only use the first two words
        exp_species = get_species_name_without_subspecies(exp_dict["species"])
        if exp_species == species:
            species_experiments.append(exp_dict)
    return species_experiments


async def parse_experiment_details(client: httpx.AsyncClient, exp_dict: dict, semaphore: asyncio.Semaphore) -> dict:
    """Parse an experiment dictionary to extract relevant information.

    Parameters
    ----------
    client : httpx.AsyncClient
        The HTTP client to use for making requests.
    exp_dict : dict
        The experiment dictionary to parse.
    semaphore : asyncio.Semaphore
        The semaphore to use for limiting the number of concurrent requests.

    Returns
    -------
    dict
        A dictionary containing the accession, description, and properties of the experiment

    """
    # getting accession and description
    accession = get_experiment_accession(exp_dict)
    description = get_experiment_description(exp_dict)
    # getting experiment metadata
    async with semaphore:
        exp_data = await get_data(client, ALL_EXP_URL + exp_dict["experimentAccession"])

    return {
        "accession": accession,
        "description": description,
        "properties": get_properties_values(exp_data),
    }


def enrich_with_searched_keywords(exp_dict: dict, keywords: list[str]) -> dict:
    """Enrich the experiment dictionary with the keywords found in the experiment.

    Parameters
    ----------
    exp_dict : dict
        The experiment dictionary to filter.
    keywords : list[str]
        The list of keywords to search for in the experiment.

    Returns
    -------
    dict | None
        The filtered experiment dictionary if keywords are found, otherwise None.

    """
    all_searchable_fields = [exp_dict["description"]] + exp_dict["properties"]
    exp_dict["found_keywords"] = keywords_in_fields(all_searchable_fields, keywords)
    if not exp_dict["found_keywords"]:
        exp_dict["failure_reason"] = NO_KEYWORD_FOUND_FAILURE_REASON
    return exp_dict


def extract_metadata_of_interest(platform_experiments: list[dict], experiment_details: list[dict]) -> list[dict]:
    """Get the metadata of interest for the list of experiments.

    Parameters
    ----------
    platform_experiments : list[dict]
        The list of experiments from the platform.
    experiment_details : list[dict]
        The list of experiment details.

    Returns
    -------
    list[dict]
        The filtered list of experiments.

    """
    metadata = []
    accession_to_metadata_of_interest = {
        exp_dict["experimentAccession"]: {renamed_key: exp_dict.get(original_key) for original_key, renamed_key in METADATA_KEYS.items()}
        for exp_dict in platform_experiments
    }
    for exp_dict in experiment_details:
        accession = exp_dict["accession"]
        if accession not in accession_to_metadata_of_interest:
            logger.warning(f"Accession {accession} not found in experiments.")
            continue
        metadata.append({**exp_dict, **accession_to_metadata_of_interest[accession]})
    return metadata


def sample_experiments_randomly(
    experiments: list[dict],
    sampling_size: int,
    seed: int,
) -> tuple[list[str], bool]:
    """Sample experiments randomly.

    Parameters
    ----------
    experiments : list[dict]
        The list of experiments to sample.
    sampling_size : int
        The number of samples to sample.
    seed : int
        The seed to use for the random sampling.

    Returns
    -------
    tuple[list[str], bool]
        The list of sampled experiments and a boolean indicating if the sampling quota was reached.

    """
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
    """Format species name for Expression Atlas."""
    return species.replace("_", " ").capitalize().strip()


def write_dicts_to_csv(
    list_of_dicts: list[dict],
    output_file: str,
    separator: str = ",",
    write_header: bool = True,
) -> None:
    """Write a list of dictionaries to a delimited text file.

    Parameters
    ----------
    list_of_dicts : list[dict]
        The rows to write. Dictionaries may have heterogeneous keys:
        missing keys are written as empty cells, and keys absent from
        the first row are appended as extra columns.
    output_file : Path
        The output file path.
    separator : str
        The cell separator
    write_header : bool
        Whether to write a header line with the column names.

    Raises
    ------
    ValueError
        If the list is empty, or the separator is not a single character.

    """
    if not list_of_dicts:
        raise ValueError("Cannot write an empty list: no columns.")

    # column names: keys of the first row, plus any key appearing only
    # in later rows (heterogeneous dicts get a full set of columns)
    fieldnames = list(list_of_dicts[0].keys())
    for row in list_of_dicts[1:]:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with open(output_file, "w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(
            fout,
            fieldnames=fieldnames,
            delimiter=separator,
            restval="",                     # missing keys -> empty cell
            extrasaction="ignore",          # safety net, shouldn't trigger given the union above
        )
        if write_header:
            writer.writeheader()
        writer.writerows(list_of_dicts)


# -------------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------------


def main():
    args = parse_args()

    keywords = [keyword.strip() for keyword in args.keywords.split(",") if keyword.strip()]

    selected_accessions = []
    selected_experiments = []

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSING EXPRESSION ATLAS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Getting arguments
    species_name = format_species_name(args.species)

    logger.info("Getting all Expression Atlas experiments")
    all_experiments = get_eatlas_experiments()
    logger.info(f"Found {len(all_experiments)} experiments")

    logger.info("Filtering on species name")
    species_experiments = get_species_experiments(all_experiments, species_name)
    logger.info(f"Found {len(species_experiments)} experiments for species {species_name}")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FOR DEV PURPOSES / TESTING: RESTRICT TO SPECIFIC ACCESSIONS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.accessions:
        logger.info(f"Keeping only accessions {args.accessions}")
        dev_accessions = args.accessions.split(",")
        species_experiments = [d for d in species_experiments if d["experimentAccession"] in dev_accessions]
        logger.info(f"Kept {len(species_experiments)} datasets for dev / testing purposes")

    logger.info("Filtering experiments based on platform")
    platform_experiments = filter_by_platform(species_experiments, args.platform)
    logger.info(f"{len(platform_experiments)} experiments passed filter")

    logger.info("Fetching and parsing metadata of experiments")
    experiment_details = asyncio.run(parse_details_of_all_experiments(platform_experiments))

    experiments = extract_metadata_of_interest(platform_experiments, experiment_details)

    # getting accessions of selected experiments so far
    selected_accessions = [exp_dict["accession"] for exp_dict in experiments]

    if keywords:
        logger.info(f"Filtering experiments with keywords {keywords}")
        experiments = [enrich_with_searched_keywords(exp_dict, keywords) for exp_dict in experiments]
        selected_accessions = [exp_dict["accession"] for exp_dict in experiments if exp_dict["found_keywords"]]
        logger.info(f"Found {len(selected_accessions)} experiments corresponding to keywords {keywords}")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # SAMPLING DATASETS IF NECESSARY
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    sampling_status = "ok"
    if selected_accessions and args.random_sampling_size is not None and args.random_sampling_seed is not None:
        selected_accession_to_nb_samples = [
            {
                "accession": exp_dict["accession"],
                "nb_samples": exp_dict["mumber_of_samples"],
            }
            for exp_dict in experiments
            if exp_dict["accession"] in selected_accessions
        ]

        logger.info("Sampling experiments randomly")
        selected_accessions, sampling_quota_reached = sample_experiments_randomly(
            selected_accession_to_nb_samples,
            args.random_sampling_size,
            args.random_sampling_seed,
        )
        logger.info(f"Kept {len(selected_accessions)} experiments after random sampling")
        for exp_dict in experiments:
            if "failure_reason" not in exp_dict and exp_dict["accession"] not in selected_accessions:
                exp_dict["failure_reason"] = REJECTED_FOR_SUBSAMPLING_FAILURE_REASON

        if sampling_quota_reached:
            sampling_status = "full"

    # writing status to file
    # so that the wrapper module can get the status
    with open(SAMPLING_QUOTA_OUTFILE, "w") as fout:
        fout.write(sampling_status)

    if not selected_accessions:
        logger.warning(f"Could not find experiments for species {species_name} and keywords {keywords}")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXPORTING DATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # getting metadata of selected accessions
    selected_experiments = [exp_dict for exp_dict in experiments if exp_dict["accession"] in selected_accessions]
    rejected_experiments = [exp_dict for exp_dict in experiments if exp_dict["accession"] not in selected_accessions]

    # exporting list of accessions
    logger.info(f"Writing accessions to {ACCESSION_OUTFILE_NAME}")
    with open(ACCESSION_OUTFILE_NAME, "w") as fout:
        fout.writelines([f"{acc}\n" for acc in selected_accessions])

    if selected_experiments:
        logger.info(
            f"Writing metadata of selected experiments to {SELECTED_EXPERIMENTS_METADATA_OUTFILE_NAME}",
        )
        write_dicts_to_csv(
            selected_experiments,
            SELECTED_EXPERIMENTS_METADATA_OUTFILE_NAME,
            separator="\t",
            write_header=True,
        )

    if rejected_experiments:
        logger.info(
            f"Writing metadata of rejected species experiments to {REJECTED_EXPERIMENTS_METADATA_OUTFILE_NAME}",
        )
        write_dicts_to_csv(
            rejected_experiments,
            REJECTED_EXPERIMENTS_METADATA_OUTFILE_NAME,
            separator="\t",
            write_header=True,
        )


if __name__ == "__main__":
    main()
