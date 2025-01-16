#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import polars as pl
from pathlib import Path
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Compute M-measure for each gene")
    parser.add_argument(
        "--file",
        type=Path,
        dest="cross_joined_file",
        required=True,
        help="File where each row contains counts for two genes",
    )
    parser.add_argument(
        "--task-attempts",
        dest="task_attempts",
        type=int,
        default=1,
        help="Number of task attempts",
    )
    return parser.parse_args()


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col
        for col in lf.collect_schema().names()
        if not col.startswith(ENSEMBL_GENE_ID_COLNAME)
    ]


def compute_ratios(file: Path, low_memory: bool) -> pl.LazyFrame:
    # getting ratios for each sample
    cross_join_lf = pl.scan_parquet(file, low_memory=low_memory)
    column_pairs = {
        col: f"{col}_other"
        for col in get_count_columns(cross_join_lf)
        if not col.endswith("_other")
    }
    return cross_join_lf.select(
        [pl.col(ENSEMBL_GENE_ID_COLNAME), pl.col(f"{ENSEMBL_GENE_ID_COLNAME}_other")]
        + [
            (pl.col(col) / pl.col(other_col)).log(base=2).alias(f"{col}_log_ratio")
            for col, other_col in column_pairs.items()
        ]
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    file = args.cross_joined_file

    low_memory = True if args.task_attempts > 1 else False
    ratios_lf = compute_ratios(file, low_memory)

    ratios_df = ratios_lf.collect()
    if len(ratios_df) == 0:
        raise ValueError(f"No output following treatment of file {str(file)}")

    outfilename = "ratios.parquet"
    ratios_df.write_parquet(outfilename)


if __name__ == "__main__":
    main()
