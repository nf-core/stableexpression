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
CANDIDATE_COUNTS_OUTFILENAME = "candidate_counts.parquet"


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
        "--nb-top-stable-genes",
        type=int,
        dest="nb_most_stable_genes",
        required=True,
        help="Number of top stable genes to show",
    )
    parser.add_argument(
        "--min-pct-quantile-expr-level",
        type=float,
        dest="min_pct_quantile_expr_level",
        required=True,
        help="Minimum percentage of quantile expression level",
    )
    return parser.parse_args()


def parse_stats(file: Path) -> pl.DataFrame:
    return pl.read_csv(file).select(
        pl.col(config.GENE_ID_COLNAME).cast(pl.String()),
        pl.exclude(config.GENE_ID_COLNAME).cast(pl.Float64()),
    )


def get_best_candidates(
    stat_df: pl.DataFrame,
    candidate_selection_descriptor: str,
    nb_most_stable_genes: int,
) -> list[str]:
    logger.info("Getting best candidates")
    column_for_sorting = config.SCORING_BASE_TO_STABILITY_SCORE_COLUMN[
        candidate_selection_descriptor
    ]
    return (
        stat_df.sort(column_for_sorting, descending=False, nulls_last=True)
        .head(nb_most_stable_genes)
        .select(config.GENE_ID_COLNAME)
        .to_series()
        .to_list()
    )


"""
def filter_out_genes_with_zero_counts(stat_lf: pl.LazyFrame) -> pl.LazyFrame:
    # keep only genes that show no zero count (ie. count > 0 for all samples)
    return stat_lf.filter(pl.col(config.RATIO_ZEROS_COLNAME) == 0)
"""


def filter_out_low_expression_genes(
    stat_df: pl.DataFrame, min_pct_quantile_expr_level: float
) -> pl.DataFrame:
    logger.info("Filtering out low expression genes")
    max_quantile = (
        stat_df.select(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME).max().item()
    )
    return stat_df.filter(
        pl.col(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
        >= max_quantile * min_pct_quantile_expr_level
    )


def get_counts_for_candidates(file: Path, best_candidates: list[str]) -> pl.DataFrame:
    logger.info("Getting counts for candidate genes")
    return pl.read_parquet(file).filter(
        pl.col(config.GENE_ID_COLNAME).is_in(best_candidates)
    )


def export_data(filtered_count_df: pl.DataFrame):
    """Export gene expression data to CSV files."""
    logger.info(
        f"Exporting counts for candidate genes to: {CANDIDATE_COUNTS_OUTFILENAME}"
    )
    filtered_count_df.write_parquet(CANDIDATE_COUNTS_OUTFILENAME)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    stat_df = parse_stats(args.stat_file)

    # first basic filters
    stat_df = filter_out_low_expression_genes(stat_df, args.min_pct_quantile_expr_level)
    # stat_lf = filter_out_genes_with_zero_counts(stat_lf)

    # get base candidate genes based on the chosen statistical descriptor (cv, rcvm)
    best_candidates = get_best_candidates(
        stat_df, args.candidate_selection_descriptor, args.nb_most_stable_genes
    )

    # get counts for candidate genes
    candidate_gene_count_lf = get_counts_for_candidates(
        args.count_file, best_candidates
    )

    export_data(candidate_gene_count_lf)


if __name__ == "__main__":
    main()
