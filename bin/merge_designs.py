#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILENAME = "whole_design.csv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge count datasets"
    )
    parser.add_argument(
        "--designs", type=str, dest="design_files", required=True, help="Design files"
    )
    return parser.parse_args()


def merge_designs(design_files):
    dfs = [pd.read_csv(file) for file in design_files]
    return pd.concat(dfs, ignore_index=True)


#####################################################
# EXPORT
#####################################################


def export_data(design_df: pd.DataFrame ):
    logger.info(f"Exporting normalised counts to: {OUTFILENAME}")
    design_df.to_csv(OUTFILENAME, index=False, header=True)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    design_files = [Path(file) for file in args.design_files.split(" ")]

    # putting all designs into a single dataframe
    design_df = merge_designs(design_files)
    export_data(design_df)


if __name__ == "__main__":
    main()
