#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import polars as pl
from pathlib import Path
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
RATIOS_STDS_COLNAME = "ratios_stds"

# experimentally chosen
RATIO_CHUNK_SIZE = 100


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
        dest="ratio_file",
        required=True,
        help="File log of pairwise expression ratios",
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


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col
        for col in lf.collect_schema().names()
        if not col.startswith(ENSEMBL_GENE_ID_COLNAME)
    ]


def compute_standard_deviations(file: Path, low_memory: bool) -> pl.LazyFrame:
    ratios_lf = pl.scan_parquet(file, low_memory=low_memory)
    ratio_columns = [
        col for col in ratios_lf.collect_schema().names() if col.endswith("_log_ratio")
    ]
    concat_ratios_lf = ratios_lf.select(
        [
            pl.concat_list(
                [pl.col(col) for col in ratio_columns[i : i + RATIO_CHUNK_SIZE]]
            ).alias(f"concat_list_chunk_{i // RATIO_CHUNK_SIZE}")
            for i in range(0, len(ratio_columns), RATIO_CHUNK_SIZE)
        ]
    ).select(pl.concat_list(pl.all()).alias("ratios"))
    return pl.concat(
        [
            concat_ratios_lf.select("ratios"),
            ratios_lf.select(
                pl.exclude("^.*_log_ratio$")
            ),  # ensembl_gene_id & ensembl_gene_id_other
        ],
        how="horizontal",
    ).select(
        pl.col("ratios").list.std(ddof=0).alias(RATIOS_STDS_COLNAME),
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        pl.col(f"{ENSEMBL_GENE_ID_COLNAME}_other"),
    )


def get_column_standard_deviations(std_lf: pl.LazyFrame, column: str) -> pl.LazyFrame:
    # column is either ENSEMBL_GENE_ID_COLNAME or f"{ENSEMBL_GENE_ID_COLNAME}_other"
    return (
        std_lf.group_by(column)
        .agg(RATIOS_STDS_COLNAME)  # getting list of ratio std for this gene
        .select(
            pl.col(column).alias(ENSEMBL_GENE_ID_COLNAME), pl.col(RATIOS_STDS_COLNAME)
        )
    )


def group_standard_deviations(std_lf: pl.LazyFrame) -> pl.LazyFrame:
    # getting the standard devs for genes in the ensembl_gene_id column
    std_a = get_column_standard_deviations(std_lf, column=ENSEMBL_GENE_ID_COLNAME)
    # getting the standard devs for genes in the ensembl_gene_id_other column
    std_b = get_column_standard_deviations(
        std_lf, column=f"{ENSEMBL_GENE_ID_COLNAME}_other"
    )
    # concatenating both dataframes vertically
    # if both lists of gene ids are the identical,
    # we need to collect values only for one column to avoid duplicates
    return pl.concat([std_a, std_b], how="vertical").unique(
        subset=ENSEMBL_GENE_ID_COLNAME
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    file = args.ratio_file

    low_memory = True if args.task_attempts > 1 else False
    std_lf = compute_standard_deviations(file, low_memory)
    std_lf = group_standard_deviations(std_lf)

    # when the ratio file corresponds to the same gene ids cross joined with themselves (i == i)
    # then we want only only one row per gene id

    std_df = std_lf.collect()
    if len(std_df) == 0:
        raise ValueError(f"No output following treatment of file {str(file)}")

    outfile = args.ratio_file.name.replace("ratios", "std")
    std_df.write_parquet(outfile)


if __name__ == "__main__":
    main()
