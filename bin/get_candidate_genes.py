#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
CANDIDATE_COUNTS_OUTFILENAME = "section_{}.candidate_counts.parquet"


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
        help="File containing counts for all genes",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        dest="stat_file",
        required=True,
        help="File containing statistics of expression over all datasets",
    )
    parser.add_argument(
        "--candidate_selection_descriptor",
        type=str,
        dest="candidate_selection_descriptor",
        required=True,
        help="Statistical descriptor for gene candidate selection.",
    )
    parser.add_argument(
        "--nb-candidates-per-section",
        type=int,
        dest="nb_candidates_per_section",
        required=True,
        help="Number of candidates per section to select for subsequent steps",
    )
    parser.add_argument(
        "--nb-sections",
        type=int,
        dest="nb_sections",
        required=True,
        help="Number of sections to divide the data into",
    )
    return parser.parse_args()


def parse_stats(file: Path) -> pl.DataFrame:
    return pl.read_csv(file).select(
        pl.col(config.GENE_ID_COLNAME).cast(pl.String()),
        pl.exclude(config.GENE_ID_COLNAME).cast(pl.Float64()),
    )


def add_sections(stat_df: pl.DataFrame, col: str, nb_sections: int):
    """
    Compute the quantile intervals relatively to col.
    The function assigns to each gene a quantile interval.
    """
    return stat_df.with_columns(
        (
            pl.col(col).rank(method="ordinal") / pl.col(col).count() * nb_sections
            + pl.lit(1)
        )
        .floor()
        .cast(pl.Int8)
        # we want the only value at nb_sections +1 to be nb_sections
        .replace({nb_sections + 1: nb_sections})
        .alias("section")
    ).sort(col, descending=False, nulls_last=True)


def get_best_candidates(
    stat_df: pl.DataFrame, nb_candidates_per_section: int
) -> pl.DataFrame:
    return stat_df.group_by("section", maintain_order=True).agg(
        pl.col(config.GENE_ID_COLNAME).head(nb_candidates_per_section)
    )


def get_counts_for_candidates(file: Path, best_candidates: list[str]) -> pl.DataFrame:
    return pl.read_parquet(file).filter(
        pl.col(config.GENE_ID_COLNAME).is_in(best_candidates)
    )


def export_data(df: pl.DataFrame, section: int):
    outfile = CANDIDATE_COUNTS_OUTFILENAME.format(section)
    df.write_parquet(outfile)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    stat_df = parse_stats(args.stat_file)

    # first basic filters
    # stat_df = filter_out_low_expression_genes(stat_df, args.min_pct_quantile_expr_level)
    # stat_lf = filter_out_genes_with_zero_counts(stat_lf)

    column_for_sorting = config.SCORING_BASE_TO_STABILITY_SCORE_COLUMN[
        args.candidate_selection_descriptor
    ]

    logger.info("Getting sections")
    stat_df = add_sections(stat_df, column_for_sorting, args.nb_sections)

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
        candidate_gene_count_lf = get_counts_for_candidates(
            args.count_file, best_candidates
        )
        export_data(candidate_gene_count_lf, section)


if __name__ == "__main__":
    main()
