#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
make_cross_join.py - Make cross join between 2 count dataframes
"""


import argparse
import logging
from pathlib import Path

import polars as pl
from common import export_parquet, get_nb_rows


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--file1",
        type=Path,
        dest="count_file_1",
        required=True,
        help="Chunk count file 1",
    )
    parser.add_argument(
        "--file2",
        type=Path,
        dest="count_file_2",
        required=True,
        help="Chunk count file 2",
    )
    parser.add_argument(
        "--index1",
        type=Path,
        dest="count_file_1_index",
        required=True,
        help="Index of chunk count file 1",
    )
    parser.add_argument(
        "--index2",
        type=Path,
        dest="count_file_2_index",
        required=True,
        help="Index of chunk count file 2",
    )
    return parser.parse_args()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    lf = pl.scan_parquet(args.count_file_1)
    lf_other = pl.scan_parquet(args.count_file_2)

    logger.info("Computing cross join data")
    lf = lf.join(
        lf_other, how="cross", suffix="_other"
    )  # Perform a cross join with itself

    if len(get_nb_rows(lf)) == 0:
        raise ValueError(
            f"No output following treatment of files {args.count_file_1} and {args.count_file_2}"
        )

    outfile = f"cross_join.{args.count_file_1_index}.{args.count_file_2_index}.parquet"
    lf.sink_parquet(outfile)


if __name__ == "__main__":
    main()
