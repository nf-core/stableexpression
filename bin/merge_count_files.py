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


concat_df = pd.DataFrame()
for file in args.files:
    df = pd.read_csv(file, header=0, index_col=0)
    concat_df = concat_df.merge(df, how='outer', left_index=True, right_index=True)

concat_df.to_csv(args.outfile, index=True, header=True)
