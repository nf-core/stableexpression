#!/usr/bin/env python3

import pandas as pd
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser('Concatenate CSV files and keep one header line')
    parser.add_argument('files', type=Path, nargs='+', help='Files to concatenate')
    parser.add_argument('--outfile', type=str, required=True, help='Output filename')
    return parser.parse_args()

args = parse_args()

dfs = [pd.read_csv(file, header=0, index_col=0) for file in args.files]

concat_df = pd.concat(dfs, axis=1)

concat_df.to_csv(args.outfile, index=True, header=True)
