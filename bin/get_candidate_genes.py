#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
"""
get_candidate_genes.py - Build a shortlist of candidate genes, section by section, based on the coefficient of variation.
"""

import argparse
import logging
from pathlib import Path

import config
from common import export_parquet, parse_table

import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
CANDIDATE_COUNTS_OUTFILENAME = "section_{}.candidate_counts.parquet"
STATS_WITH_SECTION_OUTFILENAME = "section_{}.stats.parquet"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get statistics from count data for each gene"
    )
    parser.add_argument(
        "--counts",
        type=Path,
        dest="count_file",
        required=True,
        help="Files containing counts for all genes",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        dest="stat_file",
        required=True,
        help="Files containing statistics of expression over all datasets",
    )
    parser.add_argument(
        "--sections",
        type=Path,
        dest="section_file",
        required=True,
        help="File containing the section for each gene",
    )
    parser.add_argument(
        "--nb-candidates-per-section",
        type=int,
        dest="nb_candidates_per_section",
        help="Number of candidates per section to select for subsequent steps",
    )
    return parser.parse_args()


def add_sections(stat_df: pl.DataFrame, section_df: pl.DataFrame):
    """
    Add each gene's section in the dataframe containing statistics
    """
    stat_df = stat_df.join(section_df, how="left", on=config.GENE_ID_COLNAME)
    at_least_one_section_is_null = stat_df.select(pl.col(config.SECTION_COLNAME).is_null().any()).item()
    if at_least_one_section_is_null:
        raise ValueError("Section was not provided for at least one gene.")
    return stat_df


def get_best_candidates(
    stat_df: pl.DataFrame, nb_candidates_per_section: int | None
) -> pl.DataFrame:
    if nb_candidates_per_section:
        agg_expr = pl.col(config.GENE_ID_COLNAME).head(nb_candidates_per_section)
    else:
        agg_expr = pl.col(config.GENE_ID_COLNAME)
    return (
        stat_df.sort(
            config.COEFFICIENT_OF_VARIATION_COLNAME,
            descending=False,
            nulls_last=True,
            maintain_order=True,
        )
        .group_by("section", maintain_order=True)
        .agg(agg_expr)
    )


def get_counts_for_candidates(file: Path, best_candidates: list[str]) -> pl.DataFrame:
    return pl.read_parquet(file).filter(
        pl.col(config.GENE_ID_COLNAME).is_in(best_candidates)
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    stat_df = pl.read_parquet(args.stat_file)

    # adding gene sections in the stat dataframe
    logger.info("Adding sections")
    section_df = parse_table(args.section_file)
    stat_df = add_sections(stat_df, section_df)

    logger.info("Getting best candidates")
    # get base candidate genes based on the chosen statistical descriptor (cv, rcvm)
    best_candidates_df = get_best_candidates(
        stat_df,
        args.nb_candidates_per_section,
    )

    logger.info("Getting counts of best candidates")
    # this was coded as a loop in order to keep it simple
    # since it does not impact much speed and scability
    for row in best_candidates_df.iter_rows():
        section = row[0]
        best_candidates = row[1]
        count_df = get_counts_for_candidates(
            args.count_file, best_candidates
        )

        # exporting count data for the best candidates for this section
        export_parquet(
            count_df, CANDIDATE_COUNTS_OUTFILENAME.format(section)
        )

        # exporting statistics for all genes in this section
        export_parquet(
            stat_df.filter(pl.col("section") == section),
            STATS_WITH_SECTION_OUTFILENAME.format(section),
        )


if __name__ == "__main__":
    main()
