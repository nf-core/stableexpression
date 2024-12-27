#!/usr/bin/env python3

import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser("Concatenate CSV files and keep one header line")
    parser.add_argument("files", nargs="+", help="Files to concatenate")
    parser.add_argument("--outfile", type=str, required=True, help="Output filename")
    return parser.parse_args()


args = parse_args()

print("Concatenating design CSV files")
dfs = [pd.read_csv(file, header=0) for file in args.files]
concatenated_df = pd.concat(dfs, ignore_index=True)

print(f"Writing concatenated CSV file to {args.outfile}")
concatenated_df.to_csv(args.outfile, index=False)
