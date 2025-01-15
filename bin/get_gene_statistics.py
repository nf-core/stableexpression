#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
MOST_STABLE_GENES_RESULT_OUTFILENAME = "stats_most_stable_genes.csv"
ALL_GENES_RESULT_OUTFILENAME = "stats_all_genes.csv"
LOG_COUNT_SUMMARY_OUTFILENAME = "log_counts.csv"

# column names
ORIGINAL_GENE_ID_COLNAME = "original_gene_id"
ORIGINAL_GENE_IDS_COLNAME = "original_gene_ids"
ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
M_MEASURE_COLNAME = "m_measure"
GENE_NAME_COLNAME = "name"
GENE_DESCRIPTION_COLNAME = "description"
VARIATION_COEFFICIENT_COLNAME = "variation_coefficient"
STANDARD_DEVIATION_COLNAME = "standard_deviation"
MEAN_COLNAME = "mean"
EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME = "expression_level_quantile_interval"
QUANTILE_INTERVAL_STATUS_COLNAME = "quantile_interval_status"

MOST_STABLE_GENES_RESULT_COLS = [
    ENSEMBL_GENE_ID_COLNAME,
    M_MEASURE_COLNAME,
    STANDARD_DEVIATION_COLNAME,
    VARIATION_COEFFICIENT_COLNAME,
    MEAN_COLNAME,
    EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME,
    QUANTILE_INTERVAL_STATUS_COLNAME,
    GENE_NAME_COLNAME,
    GENE_DESCRIPTION_COLNAME,
    ORIGINAL_GENE_IDS_COLNAME,
]

ALL_GENES_STATS_COLS = [
    ENSEMBL_GENE_ID_COLNAME,
    M_MEASURE_COLNAME,
    MEAN_COLNAME,
    STANDARD_DEVIATION_COLNAME,
    VARIATION_COEFFICIENT_COLNAME,
]

# quantile intervals
NB_QUANTILES = 20
MAX_SELECTED_STABLE_GENES_PER_QUANTILE_INTERVAL = 10


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
        "--metadata",
        type=str,
        dest="metadata_files",
        required=True,
        help="Metadata file",
    )
    parser.add_argument(
        "--mappings", type=str, dest="mapping_files", required=True, help="Mapping file"
    )
    parser.add_argument(
        "--m-measures",
        type=Path,
        dest="m_measure_file",
        required=True,
        help="M-measure file",
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


def concat_cast_to_string_and_drop_duplicates(files: list[Path]) -> pl.LazyFrame:
    """Concatenate LazyFrames, cast all columns to String, and drop duplicates.

    The first step is to concatenate the LazyFrames. Then, the dataframe is cast
    to String to ensure that all columns have the same data type. Finally, duplicate
    rows are dropped.
    """
    lfs = get_valid_lazy_lfs(files)
    concat_lf = pl.concat(lfs)
    # dropping duplicates
    # casting all columns to String
    return concat_lf.unique().with_columns(
        [
            pl.col(column).cast(pl.String)
            for column in concat_lf.collect_schema().names()
        ]
    )


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col for col in lf.collect_schema().names() if col != ENSEMBL_GENE_ID_COLNAME
    ]


def get_counts(file: Path) -> pl.LazyFrame:
    return pl.scan_parquet(file)


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
    return concat_lf.group_by(ENSEMBL_GENE_ID_COLNAME).agg(
        pl.col(ORIGINAL_GENE_ID_COLNAME)
        .unique()
        .sort()
        .str.join(";")
        .alias(ORIGINAL_GENE_IDS_COLNAME)
    )


def preprocess_count(count_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Preprocess the count data by replacing NaN values with 0 and sorting the dataframe.
    """
    # handling NA values (genes that are not found in all datasets)
    # replace NaN values with 0
    count_lf = count_lf.fill_null(0)
    # sorting dataframe (necessary to get consistent output)
    return count_lf.sort(ENSEMBL_GENE_ID_COLNAME, descending=False)


def transform_counts_to_log_counts(count_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Transform count data to log2 scale.

    The function converts the count data to log2 scale by adding 1 to each count
    (to avoid log(0) and stabilize variance), and then taking the log base 2 of the result.
    The ENSEMBL_GENE_ID_COLNAME column is not transformed, and the transformed columns
    retain their original names.
    """
    logger.info("Getting average log2 counts")
    count_columns = get_count_columns(count_lf)
    transformed_cols = [pl.col(ENSEMBL_GENE_ID_COLNAME)] + [
        (pl.col(column) + 1).log(2).alias(column) for column in count_columns
    ]
    return count_lf.select(transformed_cols)


def get_main_statistics(count_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute log2 count descriptive statistics for each gene in the count dataframe.
    """
    logger.info("Getting descriptive statistics")
    count_columns = get_count_columns(count_lf)
    return count_lf.with_columns(
        mean=pl.concat_list(count_columns).list.mean(),
        std=pl.concat_list(count_columns).list.std(),
    ).select(
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        pl.col("mean").alias(MEAN_COLNAME),
        pl.col("std").alias(STANDARD_DEVIATION_COLNAME),
        (pl.col("std") / pl.col("mean")).alias(VARIATION_COEFFICIENT_COLNAME),
    )


def get_quantile_intervals(lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute the quantile intervals for the mean expression levels of each gene in the dataframe.

    The function assigns to each gene a quantile interval of its mean cpm compared to all genes.
    """
    logger.info("Getting average log2 cpm quantiles")
    return lf.with_columns(
        (pl.col(MEAN_COLNAME).rank() / pl.col(MEAN_COLNAME).count() * NB_QUANTILES)
        .floor()
        .cast(pl.Int8)
        # we want the only value = NB_QUANTILES to be NB_QUANTILES - 1
        # because the last quantile interval is [NB_QUANTILES - 1, NB_QUANTILES]
        .replace({NB_QUANTILES: NB_QUANTILES - 1})
        .alias(EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
    )


def compute_general_statistics(count_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute descriptive statistics and quantile intervals for gene expression data.

    This function calculates the main statistics (mean, standard deviation, and
    variation coefficient) for each gene. It then computes the quantile intervals based on mean expression levels,
    assigning each gene a quantile interval of its mean counts per million (cpm) compared to all genes.
    """
    # getting expression statistics
    stat_lf = get_main_statistics(count_lf)

    # getting quantile intervals
    stat_lf = get_quantile_intervals(stat_lf)

    return stat_lf


def add_computed_statistics(
    stat_lf: pl.LazyFrame, m_measure_file: Path
) -> pl.LazyFrame:
    m_measure_lf = pl.scan_csv(m_measure_file)
    return stat_lf.join(m_measure_lf, on=ENSEMBL_GENE_ID_COLNAME, how="left")


def merge_data(
    stat_lf: pl.LazyFrame, metadata_lf: pl.LazyFrame, mapping_lf: pl.LazyFrame
) -> pl.LazyFrame:
    """Merge the statistics dataframe with the metadata dataframe and the mapping dataframe."""
    # we need to ensure that the index of stat_lf are strings
    return (
        stat_lf.join(metadata_lf, on=ENSEMBL_GENE_ID_COLNAME, how="left")
        .join(mapping_lf, on=ENSEMBL_GENE_ID_COLNAME, how="left")
        .unique()  # just in case
        .sort(STANDARD_DEVIATION_COLNAME, descending=False)  # VERY IMPORTANT
    )


def get_status(quantile_interval: int) -> str:
    """Return the expression level status of the gene given its quantile interval."""
    if quantile_interval == NB_QUANTILES - 1:
        return "very_high_expression"
    elif quantile_interval == NB_QUANTILES - 2:
        return "high_expression"
    elif quantile_interval == 0:
        return "very_low_expression"
    elif quantile_interval == 1:
        return "low_expression"
    else:
        return "ok"


def get_most_stable_genes(stat_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Extract the most stable genes from the statistics dataframe.

    This function groups the genes by their expression level quantile intervals
    and selects the top most stable genes from each quantile interval.
    The number of genes selected per interval is determined by the constant
    `MAX_SELECTED_STABLE_GENES_PER_QUANTILE_INTERVAL`.
    """
    logger.info("Getting most stable genes per quantile interval")
    return (
        stat_lf.group_by(  # the lf should be already sorted
            EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME
        )
        .head(MAX_SELECTED_STABLE_GENES_PER_QUANTILE_INTERVAL)
        .with_columns(
            pl.col(EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
            .map_elements(get_status)
            .alias(QUANTILE_INTERVAL_STATUS_COLNAME)
        )
        .sort(
            [
                EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME,
                STANDARD_DEVIATION_COLNAME,
                ENSEMBL_GENE_ID_COLNAME,  # if genes have same quantile and standard deviation, we sort them anyway by gene ID
            ],
            descending=[True, False, False],
        )
        .select(
            [
                column
                for column in MOST_STABLE_GENES_RESULT_COLS
                if column in stat_lf.collect_schema().names()
            ]
        )
    )


def format_all_genes_dataframe(stat_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Format the dataframe containing statistics for all genes by selecting the right columns
    and sorting the dataframe by gene ID.
    """
    return stat_lf.select(
        [
            column
            for column in ALL_GENES_STATS_COLS
            if column in stat_lf.collect_schema().names()
        ]
    ).sort(ENSEMBL_GENE_ID_COLNAME, descending=False)


def export_data(
    most_stable_genes_lf: pl.LazyFrame,
    all_genes_stat_lf: pl.LazyFrame,
    log_count_lf: pl.LazyFrame,
):
    """Export gene expression data to CSV files."""
    logger.info(
        f"Exporting statistics for the most stable genes to: {MOST_STABLE_GENES_RESULT_OUTFILENAME}"
    )
    most_stable_genes_lf.collect().write_csv(MOST_STABLE_GENES_RESULT_OUTFILENAME)

    logger.info(
        f"Exporting statistics for all genes to: {ALL_GENES_RESULT_OUTFILENAME}"
    )
    all_genes_stat_lf.collect().write_csv(ALL_GENES_RESULT_OUTFILENAME)

    logger.info(f"Exporting log normalised counts to: {LOG_COUNT_SUMMARY_OUTFILENAME}")
    log_count_lf.collect().write_csv(LOG_COUNT_SUMMARY_OUTFILENAME)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    metadata_files = [Path(file) for file in args.metadata_files.split(" ")]
    mapping_files = [Path(file) for file in args.mapping_files.split(" ")]

    # putting all counts into a single dataframe
    count_lf = get_counts(args.count_file)
    # preprocessing counts (removing 0 counts / log transformation)
    count_lf = preprocess_count(count_lf)

    # getting log2 counts
    log_count_lf = transform_counts_to_log_counts(count_lf)

    # getting metadata and mappings
    metadata_lf = get_metadata(metadata_files)
    mapping_lf = get_mappings(mapping_files)

    # computing statistics (mean, standard deviation, coefficient of variation, quantiles)
    stat_lf = compute_general_statistics(log_count_lf)

    # adding other statistics (example: m-measure)
    stat_lf = add_computed_statistics(stat_lf, args.m_measure_file)

    stat_lf = merge_data(stat_lf, metadata_lf, mapping_lf)

    most_stable_genes_lf = get_most_stable_genes(stat_lf)

    all_genes_stat_lf = format_all_genes_dataframe(stat_lf)
    export_data(most_stable_genes_lf, all_genes_stat_lf, log_count_lf)


if __name__ == "__main__":
    main()
