#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import polars as pl
from pathlib import Path
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
M_MEASURE_OUTFILE_NAME = "m_measures.csv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Compute M-measure for each gene")
    parser.add_argument(
        "--files",
        type=str,
        dest="std_files",
        required=True,
        help="File containing std of lof expression ratios",
    )
    parser.add_argument(
        "--task-attempts",
        dest="task_attempts",
        type=int,
        default=1,
        help="Number of task attempts",
    )
    return parser.parse_args()


def get_nb_rows(lf: pl.LazyFrame):
    return lf.select(pl.len()).collect().item()


def concat_all_std_data(files: list[Path], low_memory: bool) -> pl.LazyFrame:
    lfs = [pl.scan_parquet(file, low_memory=low_memory) for file in files]
    lf = pl.concat(lfs)
    return (
        lf.explode("ratios_std")
        .group_by(ENSEMBL_GENE_ID_COLNAME)
        .agg(pl.col("ratios_std"))
    )


def compute_m_measures(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        (
            pl.col("ratios_std").list.mean() / (pl.col("ratios_std").list.len() - 1)
        ).alias("m_measure"),
    )


##################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    files = [Path(file) for file in args.std_files.split(" ")]
    low_memory = True if args.task_attempts > 1 else False

    logger.info("Computing m-measure for all genes")
    all_std_lf = concat_all_std_data(files, low_memory)
    m_measure_lf = compute_m_measures(all_std_lf)
    m_measure_lf.collect().write_csv(M_MEASURE_OUTFILE_NAME)


if __name__ == "__main__":
    main()
