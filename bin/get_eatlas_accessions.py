#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
    before_sleep_log,
)
import json
from functools import partial
from multiprocessing import Pool
import nltk
from nltk.corpus import wordnet
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ALL_EXP_URL = "https://www.ebi.ac.uk/gxa/json/experiments/"
ACCESSION_OUTFILE_NAME = "accessions.txt"
FILTERED_EXPERIMENTS_OUTFILE_NAME = "filtered_experiments.json"

##################################################################
##################################################################
# NLTK MODELS AND OBJECTS
##################################################################
##################################################################

nltk.download("punkt_tab")
nltk.download("averaged_perceptron_tagger_eng")
nltk.download("wordnet")

lemmatizer = nltk.WordNetLemmatizer()
stemmer = nltk.PorterStemmer()

##################################################################
##################################################################
# EXCEPTIONS
##################################################################
##################################################################


class ExpressionAtlasNothingFoundError(Exception):
    pass


##################################################################
##################################################################
# FUNCTIONS
##################################################################
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Get expression atlas accessions")
    parser.add_argument("--species", type=str, help="Species to convert IDs for")
    parser.add_argument(
        "--keywords",
        type=str,
        nargs="*",
        help="Keywords to search for in experiment description",
    )
    return parser.parse_args()


def get_wordnet_pos(token: str):
    tag = nltk.pos_tag([token])[0][1][0].upper()
    tag_dict = {
        "J": wordnet.ADJ,
        "N": wordnet.NOUN,
        "V": wordnet.VERB,
        "R": wordnet.ADV,
    }
    return tag_dict.get(tag, wordnet.NOUN)  # Default to NOUN if not found


def get_stemmed_tokens(sentence: str):
    """
    Tokenize a sentence into its constituent words, and then stem each word

    Parameters
    ----------
    sentence : str
        The sentence to be tokenized and stemmed

    Returns
    -------
    tokens : List[str]
        The list of stemmed tokens
    """

    tokens = nltk.word_tokenize(sentence)
    return [stemmer.stem(token) for token in tokens]


def get_lemmed_tokens(sentence: str):
    """
    Tokenize a sentence into its constituent words, and then lemmatize each word

    Parameters
    ----------
    sentence : str
        The sentence to be tokenized and lemmatized

    Returns
    -------
    tokens : List[str]
        The list of lemmatized tokens
    """
    tokens = nltk.word_tokenize(sentence)
    return [lemmatizer.lemmatize(token, get_wordnet_pos(token)) for token in tokens]


def get_synonyms(word):
    """
    Get all synonyms of a word from the wordnet database.

    Parameters
    ----------
    word : str
        The word for which to get synonyms

    Returns
    -------
    synonyms : set
        A set of all synonyms of the word
    """
    synonyms = []
    for syn in wordnet.synsets(word):
        for lemma in syn.lemmas():
            synonyms.append(lemma.name())  # Get the name of each lemma (synonym)
    return set(synonyms)  # Return as a set to avoid duplicates


def get_all_candidate_target_words(sentence: str):
    """
    Get all candidate target words from a sentence by stemming and lemmatizing the
    tokens and getting synonyms from the wordnet database.

    Parameters
    ----------
    sentence : str
        The sentence from which to get candidate target words

    Returns
    -------
    candidates : list
        A list of all candidate target words
    """
    candidates = []
    lemmatized_tokens = get_stemmed_tokens(sentence)
    stemmed_tokens = get_stemmed_tokens(sentence)
    tokens = list(set(lemmatized_tokens + stemmed_tokens))
    for token in tokens:
        candidates += get_synonyms(token)
    return candidates


def word_in_sentence(word: str, sentence: str):
    """
    Checks if a word (or a stemmed version of it) is in a sentence, or if it is a
    subword of a stemmed version of any word in the sentence.

    Parameters
    ----------
    word : str
        The word to be searched for
    sentence : str
        The sentence in which to search for the word

    Returns
    -------
    bool
        True if the word is found in the sentence, False otherwise
    """
    for stemmed_word in [word] + get_stemmed_tokens(word):
        # testing if stemmed word is in sentence as it is
        if stemmed_word in sentence:
            return True
        # or testing if stemmed word is a subword of a stemmed word from the sentence
        for target_word in get_all_candidate_target_words(sentence):
            if stemmed_word in target_word:
                return True
    return False


@retry(
    retry=retry_if_exception_type(ExpressionAtlasNothingFoundError),
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def get_data(url: str):
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
    if response.status_code == 200:
        return response.json()
    elif response.status_code == 500:
        raise ExpressionAtlasNothingFoundError
    else:
        raise RuntimeError(f"Failed to retrieve data: {response.status_code}")


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


def get_experiment_accesssion(exp_dict: dict):
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


def get_species_experiments(
    species: str,
):
    """
    Gets all experiments for a given species

    Parameters
    ----------
    species : str
        Name of species. Example: "Arabidopsis thaliana"

    Returns
    -------
    experiments : list
        A list of experiment dictionaries
    """
    data = get_data(ALL_EXP_URL)
    experiments = []
    for exp_dict in data["experiments"]:
        if exp_dict["species"] == species:
            experiments.append(exp_dict)
    return experiments


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
    accession = get_experiment_accesssion(exp_dict)
    description = get_experiment_description(exp_dict)
    # getting properties of this experiment
    exp_data = get_experiment_data(exp_dict)
    properties_values = get_properties_values(exp_data)

    return {
        "accession": accession,
        "description": description,
        "properties": properties_values,
    }


def keywords_in_experiment(fields: list[str], keywords: list[str]):
    return [
        keyword
        for keyword in keywords
        for field in fields
        if word_in_sentence(keyword, field)
    ]


def filter_experiment(exp_dict: dict, keywords: list[str]):
    all_searchable_fields = [exp_dict["description"]] + exp_dict["properties"]
    found_keywords = keywords_in_experiment(all_searchable_fields, keywords)
    # only returning experiments if found keywords
    if found_keywords:
        exp_dict["found_keywords"] = list(set(found_keywords))
        return exp_dict
    else:
        return None


def format_species_name(species: str):
    return species.replace("_", " ").capitalize().strip()


##################################################################
##################################################################
# MAIN
##################################################################
##################################################################


def main():
    args = parse_args()

    # Getting arguments
    species_name = format_species_name(args.species)
    keywords = args.keywords

    logger.info(f"Getting experiments corresponding to species {species_name}")
    experiments = get_species_experiments(species_name)
    logger.info(f"Found {len(experiments)} experiments")

    logger.info("Parsing experiments")
    with Pool() as pool:
        results = pool.map(parse_experiment, experiments)

    if keywords:
        logger.info(f"Filtering experiments with keywords {keywords}")
        func = partial(filter_experiment, keywords=keywords)
        with Pool() as pool:
            results = [res for res in pool.map(func, results) if res is not None]

        if results:
            logger.info(f"Kept {len(results)} experiments")
        else:
            raise RuntimeError(
                f"Could not find experiments for species {args.species} and keywords {args.keywords}"
            )

    selected_accessions = [exp_dict["accession"] for exp_dict in results]
    logger.info(f"Writing accessions to {ACCESSION_OUTFILE_NAME}")
    with open(ACCESSION_OUTFILE_NAME, "w") as fout:
        fout.writelines([f"{acc}\n" for acc in selected_accessions])

    logger.info(f"Writing filtered experiments to {FILTERED_EXPERIMENTS_OUTFILE_NAME}")
    with open(FILTERED_EXPERIMENTS_OUTFILE_NAME, "w") as fout:
        json.dump(results, fout)


if __name__ == "__main__":
    main()
