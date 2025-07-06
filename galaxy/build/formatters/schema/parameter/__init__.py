from .base import BaseParameterFormatter
from .datasets import DatasetsParameterFormatter
from .required import RequiredParameterFormatter
from .default_value import DefaultValueParameterFormatter

PARAMETER_TO_CUSTOM_CLASS = {
    "datasets": DatasetsParameterFormatter,
    "normalisation_method": RequiredParameterFormatter,
    "nb_top_gene_candidates": RequiredParameterFormatter,
    "ks_pvalue_threshold": RequiredParameterFormatter,
    "species": DefaultValueParameterFormatter,
}

__all__ = ["BaseParameterFormatter"]
