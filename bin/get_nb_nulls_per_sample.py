#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

NB_NULL_VALUES_OUTFILE = "nb_null_values.csv"


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
        "--valid-gene-ids",
        type=Path,
        dest="valid_gene_ids",
        required=True,
        help="Valid gene IDs",
    )
    return parser.parse_args()


def get_counts(file: Path) -> pl.DataFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.read_parquet(file).sort(config.GENE_ID_COLNAME, descending=False)


def get_nb_valid_genes(valid_gene_ids_file: Path) -> int:
    return len(pl.read_csv(valid_gene_ids_file).to_series())


def get_nb_internal_nulls(df: pl.DataFrame) -> pl.DataFrame:
    """
    Get the number of null values per sample.
    :return:
    A polars dataframe containing 2 columns:
        - sample: name of the sample
        - nb_nulls: number of null values
    """
    return df.select(pl.exclude(config.GENE_ID_COLNAME).is_null().sum()).transpose(
        include_header=True,
        header_name=config.SAMPLE_COLNAME,
        column_names=[config.GENE_COUNT_COLNAME],
    )


def get_total_nb_nulls(df: pl.DataFrame, nb_genes: int) -> pl.DataFrame:
    return df.with_columns(
        (pl.col(config.GENE_COUNT_COLNAME) + pl.lit(nb_genes)).alias(
            config.GENE_COUNT_COLNAME
        )
    )


def export_data(df: pl.DataFrame):
    logger.info(f"Exporting statistics for all genes to: {NB_NULL_VALUES_OUTFILE}")
    df.write_csv(NB_NULL_VALUES_OUTFILE)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    # putting all counts into a single dataframe
    logger.info("Loading count data...")
    count_df = get_counts(args.count_file)
    logger.info(
        f"Loaded count data with {count_df.shape[0]} rows and {count_df.shape[1]} columns"
    )

    # getting nb of missing values inside the dataframe (rare but may exist)
    nb_null_values_df = get_nb_internal_nulls(count_df)

    # getting nb of missing valid genes inside the dataframe
    nb_valid_genes = get_nb_valid_genes(args.valid_gene_ids)
    nb_null_values_df = get_total_nb_nulls(nb_null_values_df, nb_valid_genes)

    export_data(nb_null_values_df)
    logger.info("Done")


if __name__ == "__main__":
    main()
