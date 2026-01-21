#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    file = Path(sys.argv[1])

    logger.info("Collecting statistics...")
    # parsing file manually because it's not a standard CSV format
    with open(file, "r") as f:
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

    outfile = file.name.replace(".csv", ".transposed.csv")
    logger.info(f"Saving statistics to {outfile}")
    df.T.to_csv(outfile, index=False, header=True)


if __name__ == "__main__":
    main()
