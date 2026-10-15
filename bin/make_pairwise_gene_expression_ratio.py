#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
"""
Perform a cross join of a count file with itself and compute log ratio for each pair of gene, for each column

The output is like:
    gene_id	gene_id_other	sampleA	sampleB
    AT1G01060	AT4G34290	0.48072526	-0.049063634
    AT1G01060	AT4G34460	0.21205564	0.15575644

"""

import argparse
import logging
from pathlib import Path

import config
from common import get_nb_rows, get_count_columns
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

CROSS_JOIN_INTERMEDIATE_OUTFILE = "cross_join.parquet"
JOIN_SUFFIX = "_other"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--file",
        type=Path,
        dest="count_file",
        required=True,
        help="File where each row contains counts for two genes",
    )
    return parser.parse_args()


def compute_ratios(cross_join_lf: pl.LazyFrame) -> pl.LazyFrame:
    # getting ratios for each sample

    column_pairs = {
        col: f"{col}{JOIN_SUFFIX}"
        for col in get_count_columns(cross_join_lf)
        if not col.endswith(JOIN_SUFFIX)
    }
    return cross_join_lf.select(
        [pl.col(config.GENE_ID_COLNAME), pl.col(f"{config.GENE_ID_COLNAME}{JOIN_SUFFIX}")]
        + [
            (pl.col(col) / pl.col(other_col)).log(base=2).alias(f"{col}_log_ratio")
            for col, other_col in column_pairs.items()
        ]
    )

def make_cross_join(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.join(lf, how="cross", suffix=JOIN_SUFFIX)

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    count_lf = pl.scan_parquet(args.count_file)

    logger.info(f"Computing cross join for {args.count_file}")
    cross_join_lf = make_cross_join(count_lf)

    logger.info(f"Computing ratios for {args.count_file}")
    ratios_lf = compute_ratios(cross_join_lf)

    if get_nb_rows(ratios_lf) == 0:
        raise ValueError(
            f"No output following treatment of file {args.count_file}"
        )

    outfile = f"{args.count_file.stem}.ratios.parquet"
    ratios_lf.sink_parquet(outfile)

    logger.info(f"Wrote ratios to {outfile}")


if __name__ == "__main__":
    main()
