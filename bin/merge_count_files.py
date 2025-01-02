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

    df1 = pd.read_csv(args.f1, index_col=0, header=0)
    df2 = pd.read_csv(args.f2, index_col=0, header=0)

    merged_df = df1.merge(df2, how="outer", left_index=True, right_index=True)

    merged_df.to_csv("merged.csv", index=True, header=True)
