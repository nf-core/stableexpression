#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

"""
separate_counts_by_batch.py Separate merged count dataframe into multiple dataframes, each exported dataframe corresponding to a unique batch in the design.
"""

import argparse
import logging
from pathlib import Path

import config
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# ----------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser("Rename gene IDs using mapped IDs")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Input file containing counts"
    )
    parser.add_argument(
        "--design", type=Path, dest="design_file", required=True, help="Design file"
    )
    return parser.parse_args()


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    logger.info(f"Parsing counts from {args.count_file}")
    lf = pl.scan_parquet(args.count_file)
    logger.info(f"Getting design from {args.design_file}")
    design_df = pl.read_csv(args.design_file)

    logger.info("Separating counts")
    for batch in design_df['batch'].unique():
        samples = design_df.filter(pl.col('batch') == batch)['sample'].to_list()
        batch_lf = lf.select([pl.col(config.GENE_ID_COLNAME), pl.col(samples)])
        outfile = f"{batch}.parquet"
        batch_lf.sink_parquet(outfile)

    logger.info('Done')


if __name__ == "__main__":
    main()
