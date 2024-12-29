#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import requests
from retry import retry
from os import cpu_count
import json
from multiprocessing import Pool
import nltk
from nltk.corpus import wordnet

ALL_EXP_URL = "https://www.ebi.ac.uk/gxa/json/experiments/"
ACCESSION_OUTFILE_NAME = "accessions.txt"
JSON_OUTFILE_NAME = "found.json"

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


@retry(ExpressionAtlasNothingFoundError, tries=3, delay=2, backoff=2)
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
    return values


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


def search_keywords_in_experiment(exp_dict: dict, keywords: list[str]):
    """
    Searches for keywords in an experiment's description and conditions

    Parameters
    ----------
    exp_dict : dict
        The experiment dictionary
    keywords : list[str]
        The list of keywords to search for

    Returns
    -------
    result : dict
        A dictionary with the experiment data and a description of which keyword was found
        Example: {'data': {'experiment': {...}}, 'found': {'word': 'salt', 'description': '...'}}
        If no keyword was found, returns None
    """
    exp_data = get_experiment_data(exp_dict)
    exp_description = get_experiment_description(exp_dict)

    for keyword in keywords:
        if word_in_sentence(keyword, exp_description):
            return {
                "data": exp_data,
                "found": {"word": keyword, "description": exp_description},
            }

    # if no keyword was found in the description
    # we try and find a keyword in one of the conditions of the experimental design
    exp_data = get_experiment_data(exp_dict)
    properties_values = get_properties_values(exp_data)

    properties_values_str = " ".join(properties_values)
    for keyword in keywords:
        if word_in_sentence(keyword, properties_values_str):
            return {
                "data": exp_data,
                "found": {"word": keyword, "properties": properties_values_str},
            }

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

    print(f"Getting experiments corresponding to species {species_name}")
    species_experiments = get_species_experiments(species_name)
    print(f"Found {len(species_experiments)} experiments")

    if keywords:
        print(f"Filtering experiments corresponding to keywords {keywords}")
        selected_accessions = []
        found_dict = {}
        with Pool(cpu_count()) as pool:
            items = [
                (
                    exp_dict,
                    keywords,
                )
                for exp_dict in species_experiments
            ]
            results = pool.starmap(search_keywords_in_experiment, items)
            for result in results:
                if result is not None:
                    accession = result["data"]["experiment"]["accession"]
                    selected_accessions.append(accession)
                    found_dict[accession] = result["found"]

        if not selected_accessions:
            raise RuntimeError(
                "Could not find experiments for species {args.species} and keywords {args.keywords}"
            )
        else:
            print(
                f"Kept {len(selected_accessions)} experiments:\n{selected_accessions}"
            )

        print(f"Writing logs of found keywords to {JSON_OUTFILE_NAME}")
        with open(JSON_OUTFILE_NAME, "w") as fout:
            json.dump(found_dict, fout)

    else:
        print("No keywords specified. Keeping all experiments")
        selected_accessions = [
            exp_dict["experimentAccession"] for exp_dict in species_experiments
        ]
        print(selected_accessions)

    print(f"Writing accessions to {ACCESSION_OUTFILE_NAME}")
    with open(ACCESSION_OUTFILE_NAME, "w") as fout:
        fout.writelines([f"{acc}\n" for acc in selected_accessions])


if __name__ == "__main__":
    main()
