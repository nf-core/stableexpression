import nltk
from nltk.corpus import wordnet

nltk.download("punkt_tab")
nltk.download("averaged_perceptron_tagger_eng")
nltk.download("wordnet")

lemmatizer = nltk.WordNetLemmatizer()
stemmer = nltk.PorterStemmer()


def get_wordnet_pos(token: str) -> str:
    tag = nltk.pos_tag([token])[0][1][0].upper()
    tag_dict = {
        "J": wordnet.ADJ,
        "N": wordnet.NOUN,
        "V": wordnet.VERB,
        "R": wordnet.ADV,
    }
    return tag_dict.get(tag, wordnet.NOUN)  # Default to NOUN if not found


def get_stemmed_tokens(sentence: str) -> list[str]:
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


def get_lemmed_tokens(sentence: str) -> list[str]:
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


def get_synonyms(word) -> set[str]:
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


def get_all_candidate_target_words(sentence: str) -> list[str]:
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


def word_is_in_sentence(word: str, sentence: str) -> bool:
    """
    Check if a word (or a stemmed version of it) is in a sentence, or if it is a
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


def keywords_in_fields(fields: list[str], keywords: list[str]) -> list[str]:
    return [
        keyword
        for keyword in keywords
        for field in fields
        if word_is_in_sentence(keyword, field)
    ]

