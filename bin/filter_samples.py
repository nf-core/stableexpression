#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
from common import export_parquet, parse_count_table

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".filtered.parquet"
RATIO_NULL_VALUES_PER_SAMPLE_OUTFILE = "ratio_null_values_per_sample.csv"
RATIO_ZEROS_PER_SAMPLE_OUTFILE = "ratio_zeros_per_sample.csv"
RATIO_NULL_VALUES_OUTFILE = "ratio_null_values.csv"
RATIO_ZERO_VALUES_OUTFILE = "ratio_zeros.csv"
NB_REJECTED_SAMPLES_OUTFILE = "nb_rejected_samples.csv"
NB_KEPT_SAMPLES_OUTFILE = "nb_kept_samples.csv"


#####################################################
# FUNCTIONS
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Filter out samples not valid")
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
    parser.add_argument(
        "--max-null-ratio",
        type=float,
        dest="max_null_ratio",
        required=True,
        help="Maximum ratio of null values",
    )
    parser.add_argument(
        "--max-zero-ratio",
        type=float,
        dest="max_zero_ratio",
        required=True,
        help="Maximum ratio of zeros allowed",
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
# MAIN
#####################################################


def main():
    args = parse_args()

    # putting all counts into a single dataframe
    logger.info("Loading count data...")
    count_df = parse_count_table(args.count_file)
    nb_genes = count_df.height
    nb_samples = count_df.width - 1
    logger.info(
        f"Loaded count data with {nb_genes} genes and {nb_samples} samples"
    )

     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # FILTERING OUT SAMPLES WITH TOO MANY NULLS
     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

    count_df = count_df.select([config.GENE_ID_COLNAME] + valid_samples)
    # if at least one valid sample is remaining, making an updated count dataframe
    if valid_samples:
        logger.info(f"Filtered out {count_df.width - 1 - len(valid_samples)} columns")
    else:
        logger.error("No valid columns remaining")


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FILTERING OUT SAMPLES WITH TOO MANY ZEROS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # computing the number of zeros values per sample
    ratio_zeros_df = count_df.select(
        pl.exclude(config.GENE_ID_COLNAME).eq(pl.lit(0)).mean()
    )

    # getting the samples with a zero ratio lower than the max zero ratio
    valid_samples = [
        col
        for col in ratio_zeros_df.columns
        if ratio_zeros_df[col][0] <= args.max_zero_ratio
    ]

    # if at least one valid sample is remaining, making an updated count dataframe
    if valid_samples:
        logger.info(f"Filtered out {count_df.width - 1 - len(valid_samples)} columns")
        count_df = count_df.select(
            pl.col(config.GENE_ID_COLNAME), pl.col(valid_samples)
        )

        outfilename = args.count_file.with_suffix(OUTFILE_SUFFIX).name
        export_parquet(count_df, outfilename)
    else:
        logger.error("No valid columns remaining")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # EXPORTING
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # collect all ratio values for export
    ratio_values = ratio_values_df.select(config.RATIO_COLNAME).to_series().to_list()
    with open(RATIO_NULL_VALUES_OUTFILE, "w") as outfile:
        # sorting values in order to having consistent output
        outfile.write(",".join([str(val) for val in sorted(ratio_values)]))

    ratio_values_df.write_csv(RATIO_NULL_VALUES_PER_SAMPLE_OUTFILE)

    # collect all ratio values for export
    ratio_values = list(ratio_zeros_df.row(0))
    with open(RATIO_ZERO_VALUES_OUTFILE, "w") as outfile:
        # sorting values in order to having consistent output
        outfile.write(",".join([str(val) for val in sorted(ratio_values)]))

    ratio_zeros_df.write_csv(RATIO_ZEROS_PER_SAMPLE_OUTFILE)

    with open(NB_KEPT_SAMPLES_OUTFILE, "w") as fout:
        fout.write(str(len(valid_samples)))

    with open(NB_REJECTED_SAMPLES_OUTFILE, "w") as fout:
        fout.write(str(nb_samples - len(valid_samples)))

    logger.info("Done")


if __name__ == "__main__":
    main()
