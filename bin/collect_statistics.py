#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description="Collect statistics")
    parser.add_argument(
        "--file",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--cpus", type=int, dest="nb_cpus", required=True, help="Number of CPUs"
    )
    parser.add_argument(
        "--memory", type=str, dest="memory", required=True, help="Memory in GB"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    logger.info("Collecting statistics...")
    # parsing file manually because it's not a standard CSV format
    with open(args.file, "r") as f:
        lines = f.readlines()
    data = [line.strip().split(",") for line in lines]

    # getting max number of columns
    max_nb_cols = max(len(row) for row in data)
    # fill missing values with None
    for row in data:
        row += [None] * (max_nb_cols - len(row))

    df = pd.DataFrame(data)
    # the first item is the dataset name
    df.set_index(df.columns[0], inplace=True)

    outfile = args.file.name.replace(".csv", ".transposed.csv")
    logger.info(f"Saving statistics to {outfile}")
    df.T.to_csv(outfile, index=False, header=True)


if __name__ == "__main__":
    main()
