#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

M_MEASURE_OUTFILE_NAME = "m_measures.csv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Compute M-measure for each gene")
    parser.add_argument(
        "--std-files",
        type=str,
        dest="std_files",
        required=True,
        help="File containing std of lof expression ratios",
    )
    return parser.parse_args()


def get_nb_rows(lf: pl.LazyFrame):
    return lf.select(pl.len()).collect().item()


def concat_all_std_data(files: list[Path]) -> pl.LazyFrame:
    """
    Concatenate all std data from the given files into a single LazyFrame.
    Explode the ratios_stds column to get one row per gene_id and ratio_std,
    then group by gene_id again to aggregate all ratio values per gene_id.
    Each file in files is like:
    ┌────────────────┬─────────────────────────────────┐
    │ gene_id        ┆ ratios_stds                     │
    │ ---            ┆ ---                             │
    │ str            ┆ list[f32]                       │
    ╞════════════════╪═════════════════════════════════╡
    │ PRUPE_1G033600 ┆ [14.52564, 10.279425, … 10.209… │
    │ PRUPE_1G176100 ┆ [10.240738, 10.267249, … 14.51… │
    └────────────────┴─────────────────────────────────┘
    """
    lfs = [pl.scan_parquet(file) for file in files]
    lf = pl.concat(lfs)
    return (
        lf.explode(config.RATIOS_STD_COLNAME)
        .group_by(config.GENE_ID_COLNAME)
        .agg(pl.col(config.RATIOS_STD_COLNAME))
    )


def compute_m_measures(lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute the m-measure for each gene.
    The m-measure is the sum of the ratios standard deviations divided by the number of ratios minus 1.
    """
    return lf.select(
        pl.col(config.GENE_ID_COLNAME),
        (
            pl.col(config.RATIOS_STD_COLNAME).list.sum()
            / (pl.col(config.RATIOS_STD_COLNAME).list.len() - 1)
        ).alias(config.GENORM_M_MEASURE_COLNAME),
    )


def get_chunks(lst: list, chunksize: int):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), chunksize):
        yield lst[i : i + chunksize]


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    files = [Path(file) for file in args.std_files.split(" ")]

    # parsing files and making concatenation
    concat_lf = concat_all_std_data(files)

    # sort everything
    # this is very weird, but if we do not sort the ratios_std list, the M measure computation may be slightly inconsistent
    # the ratios lists should be already sorted, but if they are not, we sort them here to ensure consistency
    concat_lf = concat_lf.sort(config.GENE_ID_COLNAME).with_columns(
        pl.col(config.RATIOS_STD_COLNAME).list.sort()
    )

    # computing M measures for these gene IDs
    m_measure_lf = compute_m_measures(concat_lf)

    if m_measure_lf.select(config.GENE_ID_COLNAME).collect().is_duplicated().any():
        raise ValueError("Duplicate values found for gene IDs!")

    m_measure_lf.sink_csv(
        M_MEASURE_OUTFILE_NAME, float_precision=config.DEFAULT_CSV_FLOAT_PRECISION
    )


if __name__ == "__main__":
    main()
