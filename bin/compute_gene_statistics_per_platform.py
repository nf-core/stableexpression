#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging

from stability_scorer import StabilityScorer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"

# outfile names
ALL_GENES_RESULT_OUTFILENAME = "stats_all_genes.csv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get base statistics from count data for each gene"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--platform", type=str, required=True, help="Platform name"
    )
    return parser.parse_args()


def get_counts(file: Path) -> pl.LazyFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.scan_parquet(file).sort(ENSEMBL_GENE_ID_COLNAME, descending=False)


def export_data(
    stat_lf: pl.LazyFrame, platform: str
):
    """Export gene expression data to CSV files."""
    outfilename = f"{platform}_{ALL_GENES_RESULT_OUTFILENAME}"
    logger.info(
        f"Exporting statistics for all genes to: {outfilename}"
    )
    stat_lf.collect().write_csv(outfilename)
    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    # putting all counts into a single dataframe
    count_lf = get_counts(args.count_file)

    # computing statistics (mean, standard deviation, coefficient of variation, quantiles)
    stability_scorer = StabilityScorer(count_lf, args.platform)
    stat_lf = stability_scorer.compute_statistics_and_score()

    # exporting computed data
    export_data(stat_lf, args.platform )


if __name__ == "__main__":
    main()
