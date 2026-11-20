"""Keyword matching over experiment metadata fields.

Three evidence tiers, by decreasing confidence:

1. Lexical rules (WordNet + MeSH + curated aliases):
   prose / vocabulary / "is-a" matches. High precision, auditable.
2. Concept linking (MeSH aliases from the descriptor file):
   covered by the lexical tier through descendant terms.
3. Semantic match (scispaCy en_core_sci_sm chunking + BioLORD-2023
   embeddings): catches paraphrases lexical rules cannot reach
   ("oncogene-induced senescence" ~ "cancer"); its matches are merged
   on top and logged separately.

Models (spaCy, BioLORD) load lazily on first semantic use, so the
lexical tier works without them.

Environment: MESH_FOLDER must point to a folder for the MeSH
descriptor file (downloaded once, ~17 MB).
"""

import csv
import gzip
import importlib.util
import logging
import os
import sys
import tarfile
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from functools import cache, lru_cache
from pathlib import Path

import httpx
import nltk
from nltk.corpus import wordnet
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.getLogger("httpx").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONSTANTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# requests
STOP_RETRY_AFTER_DELAY = 600

# Mesh
MESH_DESC_URL = "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2026.gz"
MESH_FOLDER_ENV_VAR = "MESH_FOLDER"

# Scispcay
SCISPACY_MODEL_URL = "https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_core_sci_sm-0.5.4.tar.gz"
SCISPACY_MODELS_CACHE_DIR_ENV_VAR = "SCISPACY_MODELS_CACHE_DIR"

# how many levels below a keyword's MeSH records to include in the expansion
MESH_EXPLOSION_DEPTH = None

# cosine similarity threshold for the semantic (BioLORD) match
SEMANTIC_THRESHOLD = 0.60

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STARTUP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

logger.info("Loading NLP models...")
nltk.download("punkt_tab", quiet=True)
nltk.download("averaged_perceptron_tagger_eng", quiet=True)
nltk.download("wordnet", quiet=True)

lemmatizer = nltk.WordNetLemmatizer()

# heavy models, loaded lazily by get_nlp() / get_embedder()
nlp = None
embedder = None


def get_nlp():
    """Load the scispaCy model on first use only.

    en_core_sci_sm extracts biological entities (diseases, tissues, species...)
    """
    global nlp
    if nlp is None:
        import spacy
        nlp = spacy.load("en_core_sci_sm")
    return nlp


def get_embedder():
    """Load the BioLORD embedding model on first use only.

    BioLORD-2023 is trained on biological and medical ontologies (UMLS, MeSH, taxonomy).
    """
    global embedder
    if embedder is None:
        from sentence_transformers import SentenceTransformer
        logging.info("Loading BioLORD-2023 embedding model...")
        embedder = SentenceTransformer("FremyCompany/BioLORD-2023")
        logging.info("BioLORD-2023 embedding model loaded.")
    return embedder


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# general functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_cache_folder(env_var_name: str) -> Path:
    """Get the cache folder from the environment variable."""
    cache_folder = os.getenv(env_var_name)
    if not cache_folder:
        raise RuntimeError(f"The {env_var_name} environment variable must be set.")
    return Path(cache_folder)


@retry(
    stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def download_if_not_exists(url: str, target_folder: Path) -> Path:
    """Download the MeSH descriptor file to the target path."""
    target_folder.mkdir(parents=True, exist_ok=True)
    target_file = target_folder / url.rsplit("/", maxsplit=1)[-1]
    if target_file.exists():
        return target_file

    tmp_file = target_file.with_name("tmp.gz")
    with httpx.stream("GET", url) as response:
        response.raise_for_status()
        with open(tmp_file, "wb") as fout:
            for chunk in response.iter_bytes():
                fout.write(chunk)
    # rename only after a complete download: no truncated file can survive a failure
    tmp_file.rename(target_file)
    return target_file


def extract_tarball(archive: Path, target_folder: Path) -> None:
    """Extract in a temp dir first, then move: a partial extraction is never mistaken for a complete model.

    `target_folder` must NOT exist yet: it is created by the rename.
    """
    with tempfile.TemporaryDirectory(dir=target_folder.parent) as tmp:
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(tmp)   # archive contains en_core_sci_sm-0.5.4/...
        (Path(tmp) / target_folder.name).rename(target_folder)
    archive.unlink()
    logger.info(f"{archive} extracted to {target_folder}")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# file download
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def get_mesh_descriptor_file() -> Path:
    """Return the path to the MeSH descriptor file, downloading it if needed.

    The target folder is read from the MESH_FOLDER environment variable.
    The file is downloaded only if it does not exist yet.

    Returns
    -------
    Path
        The path to the (possibly freshly downloaded) MeSH descriptor file.

    Raises
    ------
    RuntimeError
        If the MESH_FOLDER environment variable is not set.

    """
    mesh_folder = get_cache_folder(MESH_FOLDER_ENV_VAR)
    logger.info(f"Downloading MeSH descriptors from {MESH_DESC_URL}")
    mesh_desc_file = download_if_not_exists(MESH_DESC_URL, mesh_folder)
    logger.info(f"Downloaded MeSH descriptors to {mesh_desc_file}")
    return mesh_desc_file


def prepare_scispacy_model() -> None:
    """Make the en_core_sci_sm model importable, downloading it if needed.

    spaCy resolves model names through regular Python imports: if the
    package `en_core_sci_sm` is on sys.path, spacy.load("en_core_sci_sm")
    works without any pip installation. So we download the tarball,
    extract it, and insert the extraction directory into sys.path.
    """
    logger.info("Preparing scispaCy model")
    if importlib.util.find_spec("en_core_sci_sm") is not None:
        module = importlib.import_module("en_core_sci_sm")
        if hasattr(module, "load"):
            return
        logger.warning("en_core_sci_sm found on sys.path but has no load(); reinstalling from cache")

    downloads_folder = get_cache_folder(SCISPACY_MODELS_CACHE_DIR_ENV_VAR)
    model_name = SCISPACY_MODEL_URL.rsplit("/", maxsplit=1)[-1].removesuffix(".tar.gz")
    model_name_no_version, _ = model_name.rsplit("-", maxsplit=1)
    model_folder = downloads_folder / model_name

    if not (model_folder / model_name_no_version).exists():
        archive = download_if_not_exists(SCISPACY_MODEL_URL, downloads_folder)
        extract_tarball(archive, model_folder)

    sys.path.insert(0, str(model_folder))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MESH index & descendants
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def open_mesh_file(mesh_desc_file: Path):
    """Open the MeSH descriptor file for streaming, gzip or plain.

    Handle cases where the file suffix is .gz while the file is not
    compressed, and cases where it is really compressed.
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
        Maps a record ID to its (non-excluded) tree numbers. Example:
        "D009369" -> {"C04"}. Tree numbers are dot-separated position
        codes: C04 = Neoplasms, C04.588 = a direct child, ...
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

            # all terms: preferred name + entry terms (synonyms),
            # from every concept including non-preferred ones
            terms = {
                term.findtext("String").lower()
                for term in record.iter("Term")
                if term.findtext("String")
            }
            terms.add(preferred)

            # tree numbers: the record's position in the MeSH hierarchy,
            # minus the excluded branches (e.g. Cysts under Neoplasms)
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

    return dict(term_to_ids), id_to_terms, dict(id_to_trees), dict(id_to_preferred)


def _is_under_node(tree_number: str, node: str) -> bool:
    """Return whether `tree_number` is `node` itself, or strictly below it.

    The comparison is done on dot-separated components, not on raw strings,
    to avoid the string-prefix trap where "D12.776" would wrongly match
    "D12.7765.500" (a different branch sharing the same prefix).
    """
    return tree_number == node or tree_number.startswith(node + ".")


@cache
def get_mesh_descendant_terms(
    keyword: str,
    max_depth: int | None = MESH_EXPLOSION_DEPTH,
) -> frozenset[str]:
    """All MeSH terms of records that are more specific than the keyword.

    The keyword is matched against record PREFERRED names first (strict);
    only if that yields nothing do we fall back to matching against any
    entry term. This avoids ambiguous words dragging in unrelated records
    and then sweeping in their whole subtree.

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
        records, respecting the depth limit and the excluded branches.

    """
    term_to_ids, id_to_terms, id_to_trees, id_to_preferred = get_mesh_indexes()
    keyword = keyword.lower()

    # which records does the keyword refer to?
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

    # the root "nodes" are the tree numbers of these records
    nodes = {
        tree_number
        for record_id in root_ids
        for tree_number in id_to_trees.get(record_id, ())
    }
    # how deep each node sits (for the depth limit below)
    node_depths = {node: node.count(".") for node in nodes}

    # collect all records under those nodes, within the depth limit
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
# NLTK (lemmas, synonyms, hypernyms)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def get_wordnet_pos(tag: str) -> str:
    """Map a Penn Treebank POS tag to the first character lemmatize() accepts."""
    tag_dict = {
        "J": wordnet.ADJ,
        "N": wordnet.NOUN,
        "V": wordnet.VERB,
        "R": wordnet.ADV,
    }
    return tag_dict.get(tag[0].upper(), wordnet.NOUN)


@lru_cache
def get_lemmas(sentence: str) -> tuple[str, ...]:
    """Tokenize a sentence and lemmatize each token (no stemming).

    POS-tagging is done in a single batch call over all tokens, which is
    significantly faster than one nltk.pos_tag call per token.

    Returns a tuple so the result is immutable and safely cacheable.
    """
    tokens = nltk.word_tokenize(sentence)
    if not tokens:
        return tuple()
    tagged = nltk.pos_tag(tokens)
    return tuple(
        lemmatizer.lemmatize(token, get_wordnet_pos(tag))
        for token, tag in tagged
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


@lru_cache
def get_field_vocabulary(field: str) -> frozenset[str]:
    """All WordNet vocabulary reachable from a field: lemmas, base forms, synonyms and derivational roots of every token."""
    vocab = set()
    for token in get_lemmas(field):
        if not token.isalpha():  # skip punctuation and IDs like "BM-MSC-M1"
            continue
        for base in get_morphy_forms(token):
            vocab |= get_synonyms(base)
            vocab |= get_derivational_roots(base)
    return frozenset(vocab)


@lru_cache
def get_field_isa_vocabulary(field: str) -> frozenset[str]:
    """Union of the hypernym lemmas of every token of a field.

    Cached per field: this is the expensive vocabulary used by the
    'is-a' match, and it does not depend on the keyword.
    """
    isa_vocab = set()
    for token in get_lemmas(field):
        if token.isalpha():
            isa_vocab |= get_hypernym_lemmas(token)
    return frozenset(isa_vocab)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KEYWORD FORMS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


@lru_cache
def get_keyword_forms(keyword: str) -> frozenset[str]:
    """All canonical forms of a keyword.

    Four sources, unioned:
    - the keyword itself;
    - the derivational roots of its WordNet base forms
      ("leaves" -> {"leaves", "leaf", "leave"} + roots of each);
    - the MeSH descendant terms of its records, depth-limited
      ("cancer" -> {"neoplasms", "tumor", "neoplasia", ...});
    """
    forms = set()
    for base in get_morphy_forms(keyword):
        forms |= get_derivational_roots(base)
    forms.add(keyword.lower())
    forms |= get_mesh_descendant_terms(keyword)
    return frozenset(forms)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LEXICAL MATCHING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def lexical_keywords_in_fields(
    fields: list[str],
    keywords: list[str],
) -> dict[str, set[str]]:
    """Rule-based matching: prose, vocabulary, and 'is-a' matches.

    Returns
    -------
    dict[str, set[str]]
        For each matched keyword, the set of matching fields.

    """
    matched = defaultdict(set)
    for kw in keywords:

        # keyword side: the keyword itself + its expanded forms
        kw_forms = get_keyword_forms(kw)
        for field in fields:

            # direct prose match: any keyword form appears in the field text
            prose_match = [form for form in kw_forms if form in field.lower()]
            if prose_match:
                logger.info(f"Prose match ({kw} -> {field}). Forms: {prose_match}")
                matched[kw].add(field)
                continue

            # lexical match: keyword form EQUALS a field lemma or synonym (no substring)
            field_vocab = get_field_vocabulary(field)
            if kw_forms & field_vocab:
                logger.info(f"Vocabulary match ({kw} -> {field}): {set(kw_forms & field_vocab)}")
                matched[kw].add(field)
                continue

            # "is-a" match: a field word is a kind of the keyword
            # ("carcinoma" is a "cancer")
            field_isa_vocab = get_field_isa_vocabulary(field)
            if kw_forms & field_isa_vocab:
                logger.info(f"'Is-a kind of' match ({kw} -> {field}): {set(kw_forms & field_isa_vocab)}")
                matched[kw].add(field)

    return dict(matched)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SEMANTIC MATCHING (scispaCy + BioLORD)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


@lru_cache(maxsize=128)
def get_keyword_embeddings(keywords: tuple[str, ...]):
    """Encode keywords with BioLORD, preserving input order.

    Returns
    -------
    tuple[tuple[str, ...], torch.Tensor]
        The input tuple (echoed back so callers iterate over the aligned
        names) and the normalized embeddings, one row per keyword.

    """
    embeddings = get_embedder().encode(
        list(keywords),
        convert_to_tensor=True,
        normalize_embeddings=True,
        show_progress_bar=False,
    )
    return keywords, embeddings


@lru_cache(maxsize=1024)
def get_field_chunks(field: str) -> tuple[str, ...]:
    """Extract the semantically relevant segments of a field.

    We extract named entities and noun phrases rather than only the
    full text: a long title dilutes the keyword's signal and drags the
    cosine score down, while its noun phrases carry concentrated meaning.

    Cached per field: chunk extraction (spaCy parsing) is expensive and
    fields recur across keywords and calls.
    """
    doc = get_nlp()(field)
    # detected entities and noun phrases, plus the full text
    candidates = {field.strip()}
    candidates.update(
        ent.text.strip() for ent in doc.ents if len(ent.text.strip()) > 2
    )
    candidates.update(
        chunk.text.strip() for chunk in doc.noun_chunks if len(chunk.text.strip()) > 2
    )
    return tuple(candidates)


@lru_cache(maxsize=1024)
def get_chunk_embeddings(chunks: tuple[str, ...]):
    """Encode a tuple of text chunks with BioLORD, preserving order.

    Cached per chunk tuple: fields recur across calls, so their chunk
    embeddings are computed once and reused.
    """
    if not chunks:
        return None
    return get_embedder().encode(
        list(chunks),
        convert_to_tensor=True,
        normalize_embeddings=True,
        show_progress_bar=False,
    )


def semantic_keywords_in_fields(
    fields: list[str],
    keywords: list[str],
    threshold: float = SEMANTIC_THRESHOLD,
) -> dict[str, set[str]]:
    """Match keywords to fields by embedding similarity.

    Each field is split into semantically relevant chunks (named
    entities and noun phrases), each chunk is embedded with BioLORD,
    and a keyword matches a field if any of its chunks is similar
    enough (cosine similarity >= threshold) to the keyword embedding.

    Returns
    -------
    dict[str, set[str]]
        For each matched keyword, the set of matching fields.

    """
    if not fields or not keywords:
        return {}

    # 1. extract the relevant chunks (noun phrases + entities) per field
    field_chunks = []
    chunk_to_field_idx = []
    for idx, field in enumerate(fields):
        for chunk in get_field_chunks(field):
            field_chunks.append(chunk)
            chunk_to_field_idx.append(idx)

    # 2. encode keywords and chunks; keyword names are echoed back so
    #    keywords and embeddings always stay aligned
    kw_names, kw_embeddings = get_keyword_embeddings(tuple(keywords))
    chunk_embeddings = get_chunk_embeddings(tuple(field_chunks))
    if chunk_embeddings is None:
        return {}

    # 3. cosine similarity matrix (keywords x chunks)
    from sentence_transformers import util

    cosine_scores = util.cos_sim(kw_embeddings, chunk_embeddings)

    # 4. filter and map chunks back to their original field
    matched = defaultdict(set)
    for kw_idx, kw in enumerate(kw_names):
        scores = cosine_scores[kw_idx]
        valid_chunk_indices = (scores >= threshold).nonzero(as_tuple=True)[0].tolist()
        for c_idx in valid_chunk_indices:
            matched[kw].add(fields[chunk_to_field_idx[c_idx]])
            logger.debug(f"Semantic match ({kw} -> {fields[chunk_to_field_idx[c_idx]]}) score {scores[c_idx]:.3f}")

    return dict(matched)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMBINED ENTRY POINT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def keywords_in_fields(
    fields: list[str],
    keywords: list[str],
    use_semantic: bool = True,
) -> dict[str, list[str]]:
    """Match keywords to fields, combining lexical and semantic evidence.

    Lexical matches are always trusted. Semantic matches (BioLORD) are
    added on top; semantic-only matches are logged separately so they
    can be audited. Set use_semantic=False to run the lexical tier only
    (no model download, no GPU/CPU cost).

    Returns
    -------
    dict[str, list[str]]
        For each matched keyword, the sorted list of matching fields.

    """
    matched: dict[str, set[str]] = defaultdict(set)
    for kw, fields_ in lexical_keywords_in_fields(fields, keywords).items():
        matched[kw] |= set(fields_)

    if use_semantic:
        semantic_matched = semantic_keywords_in_fields(fields, keywords)
        for kw, fields_ in semantic_matched.items():
            new_fields = set(fields_) - matched.get(kw, set())
            if new_fields:
                logger.info(f"Semantic-only match ({kw}): {new_fields}")
            matched[kw] |= set(fields_)

    return {kw: sorted(v) for kw, v in matched.items()}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPORT OF DATASETS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
