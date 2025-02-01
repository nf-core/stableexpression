#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import polars as pl
from pathlib import Path
import argparse
import logging

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
    parser.add_argument(
        "--task-attempts",
        dest="task_attempts",
        type=int,
        default=1,
        help="Number of task attempts",
    )
    return parser.parse_args()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    low_memory = True if args.task_attempts > 1 else False
    lf = pl.scan_parquet(args.count_file_1, low_memory=low_memory)
    lf_other = pl.scan_parquet(args.count_file_2, low_memory=low_memory)

    logger.info("Computing cross join data")
    lf = lf.join(
        lf_other, how="cross", suffix="_other"
    )  # Perform a cross join with itself

    df = lf.collect()
    if len(df) == 0:
        raise ValueError(
            f"No output following treatment of files {str(args.count_file_1)} and {str(args.count_file_2)}"
        )

    outfile = f"cross_join.{args.count_file_1_index}.{args.count_file_2_index}.parquet"
    df.write_parquet(outfile)


if __name__ == "__main__":
    main()
