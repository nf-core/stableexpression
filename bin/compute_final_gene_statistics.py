#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging

from stability_scorer import StabilityScorer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# nb of top stable genes to select and to display at the end
DEFAULT_NB_TOP_STABLE_GENES = 1000

# outfile names
TOP_STABLE_GENE_SUMMARY_OUTFILENAME = "top_stable_genes_summary.csv"
ALL_GENES_RESULT_OUTFILENAME = "stats_all_genes.csv"
ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME = "all_counts_filtered.parquet"
TOP_STABLE_GENES_COUNTS_OUTFILENAME = "top_stable_genes_transposed_counts_filtered.csv"

# column names
RANK_COLNAME = "Rank"
ORIGINAL_GENE_ID_COLNAME = "original_gene_id"
ORIGINAL_GENE_IDS_COLNAME = "original_gene_ids"
ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
GENE_NAME_COLNAME = "name"
GENE_DESCRIPTION_COLNAME = "description"

GENE_COUNT_COLNAME = "count"
SAMPLE_COLNAME = "sample"
NB_ZEROS_COLNAME = "nb_zeros"
STABILITY_SCORE_COLNAME = "stability_score"

VARIATION_COEFFICIENT_COLNAME = "variation_coefficient"
STANDARD_DEVIATION_COLNAME = "standard_deviation"
MEAN_COLNAME = "mean"
MEDIAN_COLNAME = "median"
MAD_COLNAME = "median_absolute_deviation"
EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME = "expression_level_quantile_interval"
EXPRESSION_LEVEL_STATUS_COLNAME = "expression_level_status"
NB_NULLS_COLNAME = "total_nb_nulls"
NB_NULLS_VALID_SAMPLES_COLNAME = "nb_nulls_valid_samples"

STATISTICS_COLS = [
    RANK_COLNAME,
    ENSEMBL_GENE_ID_COLNAME,
    STABILITY_SCORE_COLNAME,
    STANDARD_DEVIATION_COLNAME,
    VARIATION_COEFFICIENT_COLNAME,
    MEAN_COLNAME,
    MEDIAN_COLNAME,
    MAD_COLNAME,
    EXPRESSION_LEVEL_STATUS_COLNAME,
    NB_NULLS_COLNAME,
    NB_NULLS_VALID_SAMPLES_COLNAME
]

# making complete list of columns to export
final_cols = []
for col in STATISTICS_COLS:
    final_cols.append(col)
    for platform in ["rnaseq", "microarray"]:
        final_cols.append(f"{platform}_{col}")
# adding gene description columns
final_cols += [GENE_NAME_COLNAME, GENE_DESCRIPTION_COLNAME, ORIGINAL_GENE_IDS_COLNAME]

ALL_GENES_STATS_COLS = [
    ENSEMBL_GENE_ID_COLNAME,
    STABILITY_SCORE_COLNAME,
    MEAN_COLNAME,
    STANDARD_DEVIATION_COLNAME,
    VARIATION_COEFFICIENT_COLNAME,
]

# quantile intervals
NB_QUANTILES = 100

NB_TOP_GENES_TO_SHOW_IN_LOG_COUNTS = 100


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
        help="Count file"
    )
    parser.add_argument(
        "--stats",
        type=str,
        dest="platform_stat_files",
        required=True,
        help="Platform stat file"
    )
    parser.add_argument(
        "--metadata",
        type=str,
        dest="metadata_files",
        required=True,
        help="Metadata file",
    )
    parser.add_argument(
        "--mappings",
        type=str,
        dest="mapping_files",
        required=True,
        help="Mapping file"
    )
    parser.add_argument(
        "--nb-top-stable-genes",
        type=int,
        dest="nb_top_stable_genes",
        required=True,
        help="Number of top stable genes to show",
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
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME)).collect_schema().names()


def cast_count_columns_to_float32(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        [pl.col(ENSEMBL_GENE_ID_COLNAME)]
        + [pl.col(column).cast(pl.Float32) for column in get_count_columns(lf)]
    )


def join_data_on_gene_id( stat_lf: pl.LazyFrame, *lfs) -> pl.LazyFrame:
    """Merge the statistics dataframe with the metadata dataframe and the mapping dataframe."""
    # we need to ensure that the index of stat_lf are strings
    for lf in lfs:
        stat_lf = stat_lf.join(lf, on=ENSEMBL_GENE_ID_COLNAME, how="left")
    return stat_lf


def get_counts(file: Path) -> pl.LazyFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.scan_parquet(file).sort(ENSEMBL_GENE_ID_COLNAME, descending=False)


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


def get_platform_statistics(platform_stat_files: list[Path]) -> pl.LazyFrame:
    """Retrieve and concatenate metadata from a list of platform-specific statistics files."""
    lf = pl.scan_csv(platform_stat_files[0])
    if len(platform_stat_files) > 1:
        for file in platform_stat_files[1:]:
            new_df = pl.scan_csv(file)
            lf = lf.join(new_df, on=ENSEMBL_GENE_ID_COLNAME, how="left")
    return lf


def sort_dataframe(lf: pl.LazyFrame) -> pl.LazyFrame:
    return (
        lf.sort(STABILITY_SCORE_COLNAME, descending=False, nulls_last=True)
        .with_row_index(name="index")
        .with_columns((pl.col("index") + 1).alias("Rank"))
        .drop("index")
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


def get_top_stable_gene_summary(
    stat_lf: pl.LazyFrame, nb_top_stable_genes: int
) -> pl.LazyFrame:
    """
    Extract the most stable genes from the statistics dataframe.
    """
    logger.info("Getting most stable genes per quantile interval")
    mapping_dict = {
        quantile_interval: get_status(quantile_interval)
        for quantile_interval in range(NB_QUANTILES)
    }
    lf = stat_lf.head(nb_top_stable_genes).with_columns(
        pl.col(EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
        .replace_strict(mapping_dict)
        .alias(EXPRESSION_LEVEL_STATUS_COLNAME)
    )

    return lf.select(
        [column for column in final_cols if column in lf.collect_schema().names()]
    )


def format_all_genes_statistics(stat_lf: pl.LazyFrame) -> pl.LazyFrame:
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
    ).sort(STABILITY_SCORE_COLNAME, descending=False)


def get_top_stable_genes_counts(
    log_count_lf: pl.LazyFrame, top_stable_genes_summary_lf: pl.LazyFrame
) -> pl.DataFrame:
    # getting list of top stable genes with their order
    top_genes_with_order = (
        top_stable_genes_summary_lf.head(NB_TOP_GENES_TO_SHOW_IN_LOG_COUNTS)
        .select(ENSEMBL_GENE_ID_COLNAME)
        .with_row_index("sort_order")
    )

    # join to get only existing genes and maintain order
    sorted_transposed_counts_df = (
        log_count_lf.join(
            top_genes_with_order,
            on=ENSEMBL_GENE_ID_COLNAME,
            how="inner"
        )
        .sort("sort_order", descending=False)
    ).collect()

    # get the actual gene names that were found (in order)
    actual_gene_names = (
        sorted_transposed_counts_df.select(ENSEMBL_GENE_ID_COLNAME)
        .to_series()
        .to_list()
    )

    return (
        sorted_transposed_counts_df
        .drop(["sort_order", ENSEMBL_GENE_ID_COLNAME])
        .transpose(column_names=actual_gene_names)
    )


def export_data(
    top_stable_genes_summary_lf: pl.LazyFrame,
    formated_stat_lf: pl.LazyFrame,
    all_counts_lf: pl.LazyFrame,
    top_stable_genes_counts_df: pl.DataFrame,
):
    """Export gene expression data to CSV files."""
    logger.info(
        f"Exporting statistics of the top stable genes to: {TOP_STABLE_GENE_SUMMARY_OUTFILENAME}"
    )
    top_stable_genes_summary_lf.collect().write_csv(TOP_STABLE_GENE_SUMMARY_OUTFILENAME)

    logger.info(
        f"Exporting statistics for all genes to: {ALL_GENES_RESULT_OUTFILENAME}"
    )
    formated_stat_lf.collect().write_csv(ALL_GENES_RESULT_OUTFILENAME)

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
    platform_stat_files = [Path(file) for file in args.platform_stat_files.split(" ")]

    count_lf = get_counts(args.count_file)

    # getting metadata and mappings
    metadata_lf = get_metadata(metadata_files)
    mapping_lf = get_mappings(mapping_files)
    platform_stat_df = get_platform_statistics(platform_stat_files)

    # computing statistics (mean, standard deviation, coefficient of variation, quantiles)
    stability_scorer = StabilityScorer(count_lf)
    stat_lf = stability_scorer.compute_statistics_and_score()

    # add gene name, description and original gene IDs
    stat_lf = join_data_on_gene_id(stat_lf, platform_stat_df, metadata_lf, mapping_lf)

    # sort genes according to the metrics present in the dataframe
    stat_lf = sort_dataframe(stat_lf)

    # getting the most stable genes
    # we don't want to exceed 1000 (for multiqc)
    nb_top_stable_genes = min(args.nb_top_stable_genes, DEFAULT_NB_TOP_STABLE_GENES)
    top_stable_genes_summary_lf = get_top_stable_gene_summary(
        stat_lf, nb_top_stable_genes
    )

    formated_stat_lf = format_all_genes_statistics(stat_lf)

    # reducing dataframe size (it is only used for plotting by MultiQC)
    count_lf = cast_count_columns_to_float32(count_lf)

    top_stable_genes_counts_df = get_top_stable_genes_counts(
        count_lf, top_stable_genes_summary_lf
    )

    # exporting computed data
    export_data(
        top_stable_genes_summary_lf,
        formated_stat_lf,
        count_lf,
        top_stable_genes_counts_df,
    )


if __name__ == "__main__":
    main()
