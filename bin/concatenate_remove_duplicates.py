#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import pandas as pd
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser("Concatenate files and remove duplicates")
    parser.add_argument("--f1", type=Path, help="Input file 1")
    parser.add_argument("--f2", type=Path, help="Input file 2")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # read and concatenate all CSV files
    input_files = [args.f1, args.f2]

    dfs = [pd.read_csv(file, header=0) for file in input_files]
    merged_df = pd.concat(dfs, ignore_index=True)

    # removes duplicates
    merged_df.drop_duplicates(inplace=True)

    merged_df.to_csv("merged.csv", index=False)
