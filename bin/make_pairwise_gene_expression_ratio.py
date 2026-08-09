#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
from common import get_nb_rows
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Compute M-measure for each gene")
    parser.add_argument(
        "--file",
        type=Path,
        dest="cross_joined_file",
        required=True,
        help="File where each row contains counts for two genes",
    )
    return parser.parse_args()


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the config.GENE_ID_COLNAME column.

    The config.GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col
        for col in lf.collect_schema().names()
        if not col.startswith(config.GENE_ID_COLNAME)
    ]


def compute_ratios(file: Path) -> pl.LazyFrame:
    # getting ratios for each sample
    cross_join_lf = pl.scan_parquet(file)
    column_pairs = {
        col: f"{col}_other"
        for col in get_count_columns(cross_join_lf)
        if not col.endswith("_other")
    }
    return cross_join_lf.select(
        [pl.col(config.GENE_ID_COLNAME), pl.col(f"{config.GENE_ID_COLNAME}_other")]
        + [
            (pl.col(col) / pl.col(other_col)).log(base=2).alias(f"{col}_log_ratio")
            for col, other_col in column_pairs.items()
        ]
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info(f"Computing ratios for {str(args.cross_joined_file)}")
    ratios_lf = compute_ratios(args.cross_joined_file)

    if get_nb_rows(ratios_lf) == 0:
        raise ValueError(
            f"No output following treatment of file {str(args.cross_joined_file)}"
        )

    outfilename = args.cross_joined_file.name.replace("cross_join", "ratios")
    ratios_lf.sink_parquet(outfilename)

    logger.info(f"Wrote ratios to {outfilename}")


if __name__ == "__main__":
    main()
