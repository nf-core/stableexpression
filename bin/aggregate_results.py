#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging

import config

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
ALL_GENE_SUMMARY_OUTFILENAME = "all_genes_summary.csv"
TOP_STABLE_GENE_SUMMARY_OUTFILENAME = "top_stable_genes_summary.csv"
ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME = "all_counts_filtered.parquet"
TOP_STABLE_GENES_COUNTS_OUTFILENAME = "top_stable_genes_transposed_counts_filtered.csv"

# nb of top stable genes to select and to display at the end
NB_TOP_STABLE_GENES = 1000
# quantile intervals
NB_QUANTILES = 100
NB_TOP_GENES_TO_SHOW_IN_BOX_PLOTS = 100

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
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--stats",
        type=Path,
        dest="stat_file",
        required=True,
        help="File containing statistics for all genes and stability scores by candidate genes",
    )
    parser.add_argument(
        "--rnaseq",
        type=Path,
        dest="rnaseq_dataset_stat_file",
        help="File containing base statistics for all genes and for all RNAseq datasets",
    )
    parser.add_argument(
        "--microarray",
        type=Path,
        dest="microarray_dataset_stat_file",
        help="File containing base statistics for all genes and for all Microarray datasets",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        dest="metadata_files",
        required=True,
        help="Metadata file",
    )
    parser.add_argument(
        "--mappings", type=str, dest="mapping_files", required=True, help="Mapping file"
    )

    return parser.parse_args()


def is_valid_lf(lf: pl.LazyFrame, file: Path) -> bool:
    """Check if a LazyFrame is valid.

    A LazyFrame is considered valid if it contains at least one row.
    """
    try:
        return not lf.limit(1).collect().is_empty()
    except FileNotFoundError:
        # strangely enough we get this error for some files existing but empty
        logger.error(f"Could not find file {str(file)}")
        return False
    except pl.exceptions.NoDataError as err:
        logger.error(f"File {str(file)} is empty: {err}")
        return False


def get_valid_lazy_lfs(files: list[Path]) -> list[pl.LazyFrame]:
    """Get a list of valid LazyFrames from a list of files.

    A LazyFrame is considered valid if it contains at least one row.
    """
    lf_dict = {file: pl.scan_csv(file) for file in files}
    return [lf for file, lf in lf_dict.items() if is_valid_lf(lf, file)]


def cast_cols_to_string(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        [pl.col(column).cast(pl.String) for column in lf.collect_schema().names()]
    )


def concat_cast_to_string_and_drop_duplicates(files: list[Path]) -> pl.LazyFrame:
    """Concatenate LazyFrames, cast all columns to String, and drop duplicates.

    The first step is to concatenate the LazyFrames. Then, the dataframe is cast
    to String to ensure that all columns have the same data type. Finally, duplicate
    rows are dropped.
    """
    lfs = get_valid_lazy_lfs(files)
    lfs = [cast_cols_to_string(lf) for lf in lfs]
    concat_lf = pl.concat(lfs)
    # dropping duplicates
    # casting all columns to String
    return concat_lf.unique()


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID column.

    The ENSEMBL_GENE_ID column contains only gene IDs.
    """
    return (
        lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME)).collect_schema().names()
    )


def cast_count_columns_to_float32(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        [pl.col(config.ENSEMBL_GENE_ID_COLNAME)]
        + [pl.col(column).cast(pl.Float32) for column in get_count_columns(lf)]
    )


def join_data_on_gene_id(stat_lf: pl.LazyFrame, *lfs: pl.LazyFrame) -> pl.LazyFrame:
    """Merge the statistics dataframe with the metadata dataframe and the mapping dataframe."""
    # we need to ensure that the index of stat_lf are strings
    for lf in lfs:
        stat_lf = stat_lf.join(lf, on=config.ENSEMBL_GENE_ID_COLNAME, how="left")
    return stat_lf


def get_counts(file: Path) -> pl.LazyFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.scan_parquet(file).sort(config.ENSEMBL_GENE_ID_COLNAME, descending=False)


def get_metadata(metadata_files: list[Path]) -> pl.LazyFrame:
    """Retrieve and concatenate metadata from a list of metadata files."""
    return concat_cast_to_string_and_drop_duplicates(metadata_files)


def get_mappings(mapping_files: list[Path]) -> pl.LazyFrame:
    concat_lf = concat_cast_to_string_and_drop_duplicates(mapping_files)
    # group by new gene IDs and gets the lis
    """Group by new gene IDs, get the list of distinct original gene IDs and convert to a string representation."""
    # t of distinct original gene IDs for each group
    # convert the list column to a string representation
    # separate the original gene IDs with a semicolon
    return concat_lf.group_by(config.ENSEMBL_GENE_ID_COLNAME).agg(
        pl.col(config.ORIGINAL_GENE_ID_COLNAME)
        .unique()
        .sort()
        .str.join(";")
        .alias(config.ORIGINAL_GENE_IDS_COLNAME)
    )


def get_status(quantile_interval: int) -> str:
    """Return the expression level status of the gene given its quantile interval."""
    if NB_QUANTILES - 5 <= quantile_interval:
        return "Very high expression"
    elif NB_QUANTILES - 10 <= quantile_interval < NB_QUANTILES - 5:
        return "High expression"
    elif 4 < quantile_interval <= 9:
        return "Low expression"
    elif quantile_interval <= 4:
        return "Very low expression"
    else:
        return "Medium range"


def add_expression_level_status(lf: pl.LazyFrame) -> pl.LazyFrame:
    logger.info("Adding expression level status")
    mapping_dict = {
        quantile_interval: get_status(quantile_interval)
        for quantile_interval in range(NB_QUANTILES)
    }
    return lf.with_columns(
        pl.col(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
        .replace_strict(mapping_dict)
        .alias(config.EXPRESSION_LEVEL_STATUS_COLNAME)
    )


def get_all_genes_summary(
    stat_summary_lf: pl.LazyFrame, *lfs: pl.LazyFrame
) -> pl.LazyFrame:
    """
    Extract the most stable genes from the statistics dataframe.
    """
    # add gene name, description and original gene IDs to statistics summary
    stat_summary_lf = join_data_on_gene_id(stat_summary_lf, *lfs)
    stat_summary_lf = add_expression_level_status(stat_summary_lf)
    return stat_summary_lf


def get_top_stable_genes_counts(
    log_count_lf: pl.LazyFrame, stat_summary_df: pl.LazyFrame
) -> pl.DataFrame:
    # getting list of top stable genes with their order
    top_genes_with_order = (
        stat_summary_df.head(NB_TOP_GENES_TO_SHOW_IN_BOX_PLOTS)
        .select(config.ENSEMBL_GENE_ID_COLNAME)
        .with_row_index("sort_order")
    )

    # join to get only existing genes and maintain order
    sorted_transposed_counts_df = (
        log_count_lf.join(
            top_genes_with_order, on=config.ENSEMBL_GENE_ID_COLNAME, how="inner"
        ).sort("sort_order", descending=False)
    ).collect()

    # get the actual gene names that were found (in order)
    actual_gene_names = (
        sorted_transposed_counts_df.select(config.ENSEMBL_GENE_ID_COLNAME)
        .to_series()
        .to_list()
    )

    return sorted_transposed_counts_df.drop(
        ["sort_order", config.ENSEMBL_GENE_ID_COLNAME]
    ).transpose(column_names=actual_gene_names)


def export_data(
    all_genes_summary_lf: pl.LazyFrame,
    top_stable_genes_summary_lf: pl.LazyFrame,
    all_counts_lf: pl.LazyFrame,
    top_stable_genes_counts_df: pl.DataFrame,
):
    """Export gene expression data to CSV files."""
    logger.info(f"Exporting statistics of all genes to: {ALL_GENE_SUMMARY_OUTFILENAME}")
    all_genes_summary_lf.collect().write_csv(ALL_GENE_SUMMARY_OUTFILENAME)

    logger.info(
        f"Exporting statistics of the top stable genes to: {TOP_STABLE_GENE_SUMMARY_OUTFILENAME}"
    )
    top_stable_genes_summary_lf.collect().write_csv(TOP_STABLE_GENE_SUMMARY_OUTFILENAME)

    logger.info(f"Exporting all counts to: {ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME}")
    all_counts_lf.collect().write_parquet(ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME)

    logger.info(
        f"Exporting counts of the top stable genes to: {TOP_STABLE_GENES_COUNTS_OUTFILENAME}"
    )
    top_stable_genes_counts_df.write_csv(TOP_STABLE_GENES_COUNTS_OUTFILENAME)

    logger.info("Done")


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    metadata_files = [Path(file) for file in args.metadata_files.split(" ")]
    mapping_files = [Path(file) for file in args.mapping_files.split(" ")]

    count_lf = get_counts(args.count_file)

    # getting data, including metadata and mappings
    all_genes_stat_summary_lf = pl.scan_csv(args.stat_file)

    platform_datasets_stat_lfs = [
        pl.scan_csv(file)
        for file in [args.rnaseq_dataset_stat_file, args.microarray_dataset_stat_file]
        if file is not None
    ]
    metadata_lf = get_metadata(metadata_files)
    mapping_lf = get_mappings(mapping_files)

    additional_data_lfs = [metadata_lf, mapping_lf] + platform_datasets_stat_lfs
    all_genes_summary_lf = get_all_genes_summary(
        all_genes_stat_summary_lf, *additional_data_lfs
    )

    top_stable_stat_summary_lf = all_genes_summary_lf.head(NB_TOP_STABLE_GENES)

    # reducing dataframe size (it is only used for plotting by MultiQC)
    count_lf = cast_count_columns_to_float32(count_lf)
    top_stable_genes_counts_df = get_top_stable_genes_counts(
        count_lf, top_stable_stat_summary_lf
    )
    # exporting computed data
    export_data(
        all_genes_summary_lf,
        top_stable_stat_summary_lf,
        count_lf,
        top_stable_genes_counts_df,
    )


if __name__ == "__main__":
    main()
