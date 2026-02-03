#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl

from resource_management import set_max_resources
from common import export_parquet, parse_count_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".nulls_filtered.parquet"
RATIO_NULL_VALUES_PER_SAMPLE_OUTFILE = "ratio_null_values_per_sample.csv"
RATIO_NULL_VALUES_OUTFILE = "ratio_null_values.csv"
NB_REJECTED_SAMPLES_OUTFILE = "nb_rejected_samples.csv"
NB_KEPT_SAMPLES_OUTFILE = "nb_kept_samples.csv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Filter out samples not valid")
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--max-null-ratio",
        type=float,
        dest="max_null_ratio",
        required=True,
        help="Maximum ratio of null values",
    )
    parser.add_argument(
        "--valid-gene-ids",
        type=Path,
        dest="valid_gene_ids",
        required=True,
        help="Valid gene IDs",
    )
    parser.add_argument(
        "--cpus", type=int, dest="nb_cpus", required=True, help="Number of CPUs"
    )
    parser.add_argument(
        "--memory", type=str, dest="memory", required=True, help="Memory in GB"
    )
    return parser.parse_args()


def get_nb_valid_genes(valid_gene_ids_file: Path) -> int:
    with open(valid_gene_ids_file, "r") as fin:
        return len(fin.readlines())


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


def get_ratio_null_values(
    df: pl.DataFrame, nb_missing_genes: int, nb_valid_genes: int
) -> pl.DataFrame:
    return df.select(
        pl.col(config.SAMPLE_COLNAME),
        (
            (pl.col(config.GENE_COUNT_COLNAME) + pl.lit(nb_missing_genes))
            / nb_valid_genes
        ).alias(config.RATIO_COLNAME),
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    set_max_resources(args.nb_cpus, args.memory, limit_polars=True)

    # putting all counts into a single dataframe
    logger.info("Loading count data...")
    count_df = parse_count_table(args.count_file)
    nb_genes = len(count_df)
    nb_samples = count_df.shape[1] - 1
    logger.info(f"Loaded count data with {nb_genes} genes and {nb_samples} samples")

    logger.info("Computing total number of nulls per sample")

    # getting nb of missing values inside the dataframe (rare but may exist)
    nb_null_values_df = get_nb_internal_nulls(count_df)

    # getting nb of missing valid genes inside the dataframe
    nb_valid_genes = get_nb_valid_genes(args.valid_gene_ids)
    nb_missing_genes = nb_valid_genes - nb_genes

    # adding the nb of missing genes to the number of null vaues for each sample
    ratio_values_df = get_ratio_null_values(
        nb_null_values_df, nb_missing_genes, nb_valid_genes
    )

    valid_samples = (
        ratio_values_df.filter(pl.col(config.RATIO_COLNAME) <= args.max_null_ratio)
        .select(pl.col(config.SAMPLE_COLNAME))
        .to_series()
        .to_list()
    )

    # if at least one valid sample is remaining, making an updated count dataframe
    if valid_samples:
        logger.info(f"Filtered out {count_df.shape[1] - len(valid_samples)} columns")
        valid_count_df = count_df.select([config.GENE_ID_COLNAME] + valid_samples)
        export_parquet(valid_count_df, args.count_file, OUTFILE_SUFFIX)
    else:
        logger.error("No valid columns remaining")

    # collect all ratio values for export
    ratio_values = ratio_values_df.select(config.RATIO_COLNAME).to_series().to_list()
    with open(RATIO_NULL_VALUES_OUTFILE, "w") as outfile:
        outfile.write(",".join([str(val) for val in ratio_values]))

    ratio_values_df.write_csv(RATIO_NULL_VALUES_PER_SAMPLE_OUTFILE)

    with open(NB_KEPT_SAMPLES_OUTFILE, "w") as fout:
        fout.write(str(len(valid_samples)))

    with open(NB_REJECTED_SAMPLES_OUTFILE, "w") as fout:
        fout.write(str(nb_samples - len(valid_samples)))

    logger.info("Done")


if __name__ == "__main__":
    main()
