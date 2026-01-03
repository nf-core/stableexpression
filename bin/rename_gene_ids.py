#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import config
import pandas as pd
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


##################################################################
# CONSTANTS
##################################################################

RENAMED_FILE_SUFFIX = ".renamed.csv"
METADATA_FILE_SUFFIX = ".metadata.csv"
MAPPING_FILE_SUFFIX = ".mapping.csv"

WARNING_REASON_FILE = "warning_reason.txt"
FAILURE_REASON_FILE = "failure_reason.txt"

UNMAPPED_FILE_SUFFIX = "unmapped.txt"
MERGED_FILE_SUFFIX = "merged.txt"
FINAL_FILE_SUFFIX = "final.txt"

##################################################################
# FUNCTIONS
##################################################################


def parse_args():
    parser = argparse.ArgumentParser("Rename gene IDs using mapped IDs")
    parser.add_argument(
        "--count-file", type=Path, required=True, help="Input file containing counts"
    )
    parser.add_argument(
        "--mappings",
        type=Path,
        required=True,
        dest="mapping_file",
        help="Mapping file containing gene IDs",
    )
    return parser.parse_args()


def parse_table(file: Path, **kwargs):
    if file.suffix == ".csv":
        return pd.read_csv(file, header=0, **kwargs)
    else:  # .tsv
        return pd.read_csv(file, header=0, sep="\t", **kwargs)


def parse_count_table(file: Path):
    # transitting to pandas dataframe helps to avoid parsing errors
    df = parse_table(file, index_col=0)
    # whatever the name of the first col, rename it to "gene_id"
    df.index.rename(config.GENE_ID_COLNAME, inplace=True)
    df.index = df.index.astype(str)
    return pl.from_pandas(df.reset_index())


##################################################################
# MAIN
##################################################################


def main():
    args = parse_args()

    logger.info(f"Converting IDs for count file {args.count_file.name}...")

    #############################################################
    # PARSING FILES
    #############################################################

    df = parse_count_table(args.count_file)

    if df.is_empty():
        msg = "COUNT FILE IS EMPTY"
        logger.warning(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)
        sys.exit(0)

    #############################################################
    # GETTING MAPPINGS
    #############################################################

    mapping_df = parse_table(args.mapping_file)
    mapping_dict = mapping_df.set_index(config.ORIGINAL_GENE_ID_COLNAME)[
        config.GENE_ID_COLNAME
    ].to_dict()

    #############################################################
    # MAPPING GENE IDS IN DATAFRAME
    #############################################################

    # IMPORTANT: KEEPING ONLY GENES THAT HAVE BEEN CONVERTED
    # filtering the DataFrame to keep only the rows where the index can be mapped
    original_nb_genes = len(df)

    rejected_df = df.filter(~pl.col(config.GENE_ID_COLNAME).is_in(mapping_dict.keys()))
    nb_unmapped_genes = len(rejected_df)

    # df = df.loc[df.index.isin(mapping_dict)]
    df = df.filter(pl.col(config.GENE_ID_COLNAME).is_in(mapping_dict.keys()))
    nb_mapped_genes = len(df)

    with open(UNMAPPED_FILE_SUFFIX, "w") as f:
        f.write(str(nb_unmapped_genes))

    if df.is_empty():
        sample_size = min(5, nb_unmapped_genes)
        example_rejected_genes = (
            rejected_df[config.GENE_ID_COLNAME].head(sample_size).to_list()
        )
        msg = f"NO GENES WERE MAPPED. EXAMPLE OF GENE IDS: {example_rejected_genes}"
        logger.error(msg)
        with open(FAILURE_REASON_FILE, "w") as f:
            f.write(msg)

        with open(MERGED_FILE_SUFFIX, "w") as f:
            f.write("0")
        with open(FINAL_FILE_SUFFIX, "w") as f:
            f.write("0")

        sys.exit(0)

    if len(df) < original_nb_genes:
        sample_size = min(5, nb_unmapped_genes)
        example_rejected_genes = (
            rejected_df[config.GENE_ID_COLNAME].head(sample_size).to_list()
        )
        msg = (
            f"{nb_mapped_genes / original_nb_genes:.2%} of genes were mapped ({nb_mapped_genes} out of {original_nb_genes}). "
            + f"Example of unmapped genes: {example_rejected_genes}"
        )
        logger.warning(msg)
        with open(WARNING_REASON_FILE, "a") as f:
            f.write(msg)
    else:
        logger.info(
            f"All genes were mapped ({nb_mapped_genes} out of {original_nb_genes})"
        )

    logger.info("Renaming gene names")
    # renaming gene names to mapped ids using mapping dict
    df = df.with_columns(
        pl.col(config.GENE_ID_COLNAME)
        .replace(mapping_dict)
        .alias(config.GENE_ID_COLNAME)
    )

    # TODO: check is there is another way to avoid duplicate gene names
    # sometimes different gene names have the same Gene ID
    # for now, we just get the mean of values, but this is not ideal

    #############################################################
    # GENE COUNT HANDLING
    #############################################################

    # handling cases where multiple genes have the same Gene ID
    # since subsequent steps in the pipeline require integer values,
    # we need to ensure that the resulting DataFrame has integer values
    logger.info("Computing mean counts for genes with duplicate IDs")
    df = df.group_by(config.GENE_ID_COLNAME, maintain_order=True).agg(
        pl.exclude(config.GENE_ID_COLNAME).mean()
    )

    nb_merged = nb_mapped_genes - len(df)
    with open(MERGED_FILE_SUFFIX, "w") as f:
        f.write(str(nb_merged))
    with open(FINAL_FILE_SUFFIX, "w") as f:
        f.write(str(len(df)))

    #############################################################
    # WRITING OUTFILES
    #############################################################
    # writing to output file

    logger.info("Writing output file")
    outfile = args.count_file.with_name(args.count_file.stem + RENAMED_FILE_SUFFIX)
    df.write_csv(outfile)

    # making dataframe for mapping (only two columns: original and new)
    mapping_df = (
        pd.DataFrame(mapping_dict, index=[0])
        .T.reset_index()  # transpose: setting keys as indexes instead of columns
        .rename(
            columns={
                "index": config.ORIGINAL_GENE_ID_COLNAME,
                0: config.GENE_ID_COLNAME,
            }
        )
    )
    mapping_file = args.count_file.with_name(args.count_file.stem + MAPPING_FILE_SUFFIX)
    mapping_df.to_csv(mapping_file, index=False, header=True)


if __name__ == "__main__":
    main()
