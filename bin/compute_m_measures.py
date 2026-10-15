#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
"""
compute_m_measure.py - Compute M-measure for each gene
"""

import argparse
import logging

import config
from common import sink_csv_with_floats
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

M_MEASURE_OUTFILE_NAME = "m_measures.csv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ratio-std-file",
        type=str,
        dest="ratio_std_file",
        required=True,
        help="File containing std of log expression ratios",
    )
    return parser.parse_args()


def compute_m_measures(lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute the m-measure for each gene.
    The m-measure is the sum of the ratios standard deviations divided by the number of ratios minus 1.
    When the number of samples is 1 (which should not happen theoretically) dividing by the number of ratios.
    """
    nb_ratios_stds = lf.select(pl.col(config.RATIOS_STD_COLNAME).list.len()).item()
    denominator = nb_ratios_stds - 1 if nb_ratios_stds > 1 else 1
    return lf.select(
        pl.col(config.GENE_ID_COLNAME),
        (pl.col(config.RATIOS_STD_COLNAME).list.sum() / denominator).alias(config.GENORM_M_MEASURE_COLNAME),
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    ratio_std_lf = pl.scan_parquet(args.ratio_std_file)

    # computing M measures for these gene IDs
    m_measure_lf = compute_m_measures(ratio_std_lf)

    if m_measure_lf.select(config.GENE_ID_COLNAME).collect().is_duplicated().any():
        raise ValueError("Duplicate values found for gene IDs!")

    sink_csv_with_floats(m_measure_lf, M_MEASURE_OUTFILE_NAME, float_precision=5)


if __name__ == "__main__":
    main()
