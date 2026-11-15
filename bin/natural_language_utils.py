import logging
import os
from pathlib import Path
import xml.etree.ElementTree as ET
import nltk
import httpx
import gzip
from nltk.corpus import wordnet
from functools import lru_cache
from collections import defaultdict
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.getLogger("httpx").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

MESH_DESC_URL = "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2026.gz"
MESH_FOLDER_ENV_VAR = "MESH_FOLDER"
STOP_RETRY_AFTER_DELAY = 600

nltk.download("punkt_tab")
nltk.download("averaged_perceptron_tagger_eng")
nltk.download("wordnet")

lemmatizer = nltk.WordNetLemmatizer()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MESH
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING)
)
def download_mesh_descriptors(mesh_desc_file: Path) -> None:
    """Download the MeSH descriptor file to the target path."""
    tmp_file = mesh_desc_file.parent / "tmp.gz"
    with httpx.stream("GET", MESH_DESC_URL) as response:
        response.raise_for_status()
        with open(tmp_file, "wb") as fout:
            for chunk in response.iter_bytes():
                fout.write(chunk)
    # renaming file upon succesful download
    tmp_file.rename(mesh_desc_file)


def get_mesh_descriptor_file() -> Path:
    """Return the path to the MeSH descriptor file, downloading it if needed.

    The target path is read from the MESH_DESC_FILE environment variable.
    The file is downloaded only if it does not exist yet.

    Returns
    -------
    Path
        The path to the (possibly freshly downloaded) MeSH descriptor file.

    Raises
    ------
    RuntimeError
        If the MESH_DESC_FILE environment variable is not set.

    """
    mesh_folder_env = os.getenv(MESH_FOLDER_ENV_VAR)
    if not mesh_folder_env:
        raise RuntimeError(
            f"The {MESH_FOLDER_ENV_VAR} environment variable must be set "
            f"to the desired folder path of the MeSH descriptor file.",
        )

    mesh_folder = Path(mesh_folder_env)
    mesh_desc_file = (mesh_folder / MESH_DESC_URL.rsplit("/", maxsplit=1)[-1])

    if mesh_desc_file.exists():
        return mesh_desc_file

    logger.info(f"Downloading MeSH descriptors to {mesh_desc_file}")
    mesh_folder.mkdir(parents=True, exist_ok=True)
    download_mesh_descriptors(mesh_desc_file)
    logger.info(f"Downloaded MeSH descriptors to {mesh_desc_file}")
    return mesh_desc_file


def open_mesh_file(mesh_desc_file: Path):
    """Open the MeSH descriptor file for streaming, gzip or plain.

    Handle cases where the file suffix is .gz while the file is not compressed,
    and cases where it is really compressed.
    """
    with open(mesh_desc_file, "rb") as probe:
        is_gzipped = probe.read(2) == b"\x1f\x8b"   # gzip magic number
    return gzip.open(mesh_desc_file, "rb") if is_gzipped else open(mesh_desc_file, "rb")


@lru_cache(maxsize=1)
def get_mesh_indexes() -> tuple[dict, dict, dict, dict]:
    """Stream the MeSH descriptor file and build all lookup tables.

    Returns
    -------
    term_to_ids : dict[str, set[str]]
        Maps every term (descriptor name or entry term, lowercase) to the
        set of record IDs it appears in. Example: "cancer" -> {"D009369"}.
    id_to_terms : dict[str, frozenset[str]]
        Maps a record ID to all its terms (lowercase). Example:
        "D009369" (Neoplasms) -> {"neoplasms", "cancer", "tumor", ...}.
    id_to_trees : dict[str, set[str]]
        Maps a record ID to its tree numbers. Example: "D009369" -> {"C04"}.
        Tree numbers are dot-separated position codes: C04 = Neoplasms,
        C04.588 = a direct child of Neoplasms, C04.588.631 = a grandchild...
    id_to_preferred : dict[str, str]
        Maps a record ID to its preferred descriptor name (lowercase).
        Example: "D009369" -> "neoplasms".

    """
    term_to_ids = defaultdict(set)
    id_to_terms = {}
    id_to_trees = defaultdict(set)
    id_to_preferred = {}

    mesh_desc_file = get_mesh_descriptor_file()
    with open_mesh_file(mesh_desc_file) as fin:
        for _, record in ET.iterparse(fin, events=("end",)):
            if record.tag != "DescriptorRecord":
                continue

            record_id = record.findtext("DescriptorUI")
            if record_id is None:
                continue

            # preferred name of the record (its own name, not a synonym)
            preferred = record.findtext("./DescriptorName/String")
            if preferred is None:
                continue
            preferred = preferred.lower()

            # all terms: preferred name + entry terms (synonyms)
            terms = {
                term.findtext("String").lower()
                for term in record.iter("Term")
                if term.findtext("String")
            }
            terms.add(preferred)

            # tree numbers: the record's position in the MeSH hierarchy
            trees = {
                tree.text
                for tree in record.iter("TreeNumber")
                if tree.text
            }

            # filling the four tables
            id_to_terms[record_id] = frozenset(terms)
            id_to_preferred[record_id] = preferred
            if trees:
                id_to_trees[record_id] = trees
            for term in terms:
                term_to_ids[term].add(record_id)

            record.clear()   # free memory as we stream

    return dict(term_to_ids), id_to_terms, dict(id_to_trees), id_to_preferred


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MESH DESCENDANTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def _is_under_node(tree_number: str, node: str) -> bool:
    """Return whether `tree_number` is under `node`.

    Return True if `tree_number` is node `node` itself, or strictly below it.

    The comparison is done on dot-separated components, not on raw strings:
    "C04.588" is under "C04", but "C04.5885.323" is NOT under "C04.5885x".
    This avoids the string-prefix trap where "D12.776" would wrongly
    match "D12.7765.500" (a different branch sharing the same prefix).
    """
    return tree_number == node or tree_number.startswith(node + ".")


@lru_cache(maxsize=None)
def get_mesh_descendant_terms(
    keyword: str,
    max_depth: int | None = 1,
) -> frozenset[str]:
    """All MeSH terms of records that are more specific than the keyword.

    The keyword is matched against record PREFERRED names first (strict);
    only if that yields nothing do we fall back to matching against any
    entry term. This avoids ambiguous words dragging in unrelated records
    (e.g. an entry-term match on a record from a completely different
    part of the hierarchy) and then sweeping in its whole subtree.

    Parameters
    ----------
    keyword : str
        The keyword to expand (e.g. "cancer").
    max_depth : int | None
        How many levels below the keyword's records to include.
        1 = direct children only (recommended starting point),
        2 = children + grandchildren, None = the whole subtree.

    Returns
    -------
    frozenset[str]
        All lowercase terms of all records at or below the keyword's
        records, respecting the depth limit. Example for ("cancer", 1):
        {"neoplasms", "cancer", "tumor", "neoplasia", "cysts", ...}
        — the major categories of cancers, not the thousands of
        organ-specific / histologic-type leaves deeper down.

    """
    term_to_ids, id_to_terms, id_to_trees, id_to_preferred = get_mesh_indexes()
    keyword = keyword.lower()

    # strict pass: records whose PREFERRED name is the keyword
    root_ids = {
        record_id
        for record_id, preferred in id_to_preferred.items()
        if preferred == keyword
    }

    # fallback pass: records where the keyword is only an entry term
    if not root_ids:
        root_ids = term_to_ids.get(keyword, set())
    if not root_ids:
        return frozenset()

    # ── step 2: where do these records sit in the MeSH tree? ──
    # the root "nodes" are their tree numbers
    nodes = {
        tree_number
        for record_id in root_ids
        for tree_number in id_to_trees.get(record_id, ())
    }
    # how deep each node sits (for the depth limit below)
    node_depths = {node: node.count(".") for node in nodes}

    # ── step 3: collect all records under those nodes, within the depth limit ──
    matched_terms = {keyword}   # the keyword itself always matches its own terms
    for record_id, record_trees in id_to_trees.items():
        for tree_number in record_trees:
            # does this record sit under one of our nodes, dot-wise?
            for node in nodes:
                if not _is_under_node(tree_number, node):
                    continue
                # is it close enough to the surface?
                if max_depth is not None:
                    depth = tree_number.count(".") - node_depths[node]
                    if depth > max_depth:
                        continue
                matched_terms |= id_to_terms[record_id]
                break   # one qualifying tree number is enough per record
            else:
                continue
            break   # ...and one qualifying record path is enough too

    return frozenset(matched_terms)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NLTK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def get_wordnet_pos(token: str) -> str:
    """Map POS tag to first character lemmatize() accepts."""
    tag = nltk.pos_tag([token])[0][1][0].upper()
    tag_dict = {
        "J": wordnet.ADJ,
        "N": wordnet.NOUN,
        "V": wordnet.VERB,
        "R": wordnet.ADV,
    }
    return tag_dict.get(tag, wordnet.NOUN)


@lru_cache
def get_lemmas(sentence: str) -> tuple[str, ...]:
    """Tokenize a sentence and lemmatize each token (no stemming).

    Returns a tuple so the result is immutable and safely cacheable.
    """
    return tuple(
        lemmatizer.lemmatize(token, get_wordnet_pos(token))
        for token in nltk.word_tokenize(sentence)
    )


@lru_cache
def get_morphy_forms(word: str) -> frozenset[str]:
    """Get all WordNet base forms of a word, for every part of speech.

    Example: "leaves" -> {"leaf", "leave"}
    """
    return frozenset(
        base
        for pos in ("n", "v", "a", "r", "s")
        if (base := wordnet.morphy(word.lower(), pos)) is not None
    ) | {word.lower()}


@lru_cache
def get_synonyms(word: str) -> frozenset[str]:
    """Get all lemma names of the synsets of a word."""
    return frozenset(
        lemma.name()
        for synset in wordnet.synsets(word.lower())
        for lemma in synset.lemmas()
    )


@lru_cache
def get_derivational_roots(word: str) -> frozenset[str]:
    """Get derivationally related forms of a word.

    This is the bridge between word families:
    "leafy" -> {"leaf"}, "ovarian" -> {"ovary"}, "mitotic" -> {"mitosis"}
    """
    roots = set()
    for synset in wordnet.synsets(word.lower()):
        for lemma in synset.lemmas():
            roots.add(lemma.name())
            for related in lemma.derivationally_related_forms():
                roots.add(related.name())
    return frozenset(roots)


@lru_cache
def get_hypernym_lemmas(word: str) -> frozenset[str]:
    """All lemma names of all hypernyms (transitively) of a word.

    Example: "carcinoma" -> {"cancer", "malignant neoplastic disease", "tumor", ...}
    """
    lemmas = set()
    for synset in wordnet.synsets(word.lower()):
        for hypernym in synset.closure(lambda s: s.hypernyms()):
            lemmas |= {lemma.name() for lemma in hypernym.lemmas()}
    return frozenset(lemmas)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KEYWORD MATCHING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


@lru_cache
def get_keyword_forms(keyword: str) -> frozenset[str]:
    """All canonical forms of a keyword: itself, its WordNet base forms, and the derivational roots of those.

    Example: "leaves" -> {"leaves", "leaf", "leave"}  (+ roots of each)
    """
    forms = set()
    for base in get_morphy_forms(keyword):
        forms |= get_derivational_roots(base)
    forms.add(keyword.lower())
    forms |= get_mesh_descendant_terms(keyword)
    return frozenset(forms)


def keywords_in_fields(fields: list[str], keywords: list[str]) -> dict[str, list[str]]:
    """Return a dict of keywords and the fields that contain them."""
    matched = defaultdict(set)
    for kw in keywords:

        # keyword side: the keyword itself + its WordNet lemmas (leaves -> leaf, leave)
        kw_forms = get_keyword_forms(kw)
        for field in fields:

            # "is-a" match: a field word is a kind of the keyword ("carcinoma" is a "cancer")
            field_isa_vocab = set()
            for token in get_lemmas(field):
                if token.isalpha():
                    field_isa_vocab |= get_hypernym_lemmas(token)
            if field_isa_vocab & kw_forms:
                logger.info(f"'Is-a kind of' match ({kw} -> {field}) ({field_isa_vocab & kw_forms})")
                matched[kw].add(field)

    return {kw: sorted(v) for kw, v in matched.items()}
