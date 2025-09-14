#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import sys
import polars as pl
from pathlib import Path
import logging

import config

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME = "cleaned_counts_filtered.parquet"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Clean data by removing aberrant samples and performing some other cleaning operations."
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--ks-stats",
        type=Path,
        dest="ks_stats_file",
        required=True,
        help="KS stats file",
    )
    parser.add_argument(
        "--ks-pvalue-threshold",
        type=str,
        dest="ks_pvalue_threshold",
        required=True,
        help="KS p-value threshold",
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
    """Get all column names except the config.ENSEMBL_GENE_ID_COLNAME column.

    The config.ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME)).collect_schema().names()


def get_counts(
    file: Path,
) -> pl.LazyFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.scan_parquet(file).sort(config.ENSEMBL_GENE_ID_COLNAME, descending=False)


def remove_samples_with_low_ks_pvalue(
    count_lf: pl.LazyFrame, ks_stats_file: Path, ks_pvalue_threshold: str
) -> pl.LazyFrame:

    ks_stats_df = pl.read_csv(ks_stats_file, has_header=True).select([config.SAMPLE_COLNAME, config.KS_TEST_COLNAME])

    # parsing threshold
    try:
        ks_pvalue_threshold = float(ks_pvalue_threshold)
    except ValueError:
        raise ValueError(
            f"KS p-value threshold {ks_pvalue_threshold} could not be cast to float"
        )

    # logging number of samples excluded from analysis
    not_valid_samples = ks_stats_df.filter(
        ks_stats_df[config.KS_TEST_COLNAME] <= ks_pvalue_threshold
    )[config.SAMPLE_COLNAME].to_list()

    if not_valid_samples:
        logger.warning(
            f"Excluded {len(not_valid_samples)} samples showing a KS p-value below {ks_pvalue_threshold}"
        )
    else:
        logger.info("No sample was excluded")

    # getting samples for which the Kolmogorov-Smirnov test pvalue is above the threshold
    valid_samples = ks_stats_df.filter(
        ks_stats_df[config.KS_TEST_COLNAME] > ks_pvalue_threshold
    )[config.SAMPLE_COLNAME].to_list()

    if not valid_samples:
        logger.error("No more valid sample to process...")
        sys.exit(101)

    # filtering the count dataframe to keep only the valid samples
    return count_lf.select([config.ENSEMBL_GENE_ID_COLNAME] + valid_samples)


def export_data( all_counts_lf: pl.LazyFrame):
    all_counts_lf.collect().write_parquet(ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME)
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

    # removing aberrant samples (ks p-value under the threshold)
    count_lf = remove_samples_with_low_ks_pvalue(count_lf, args.ks_stats_file, args.ks_pvalue_threshold)

    # exporting computed data
    export_data(count_lf)


if __name__ == "__main__":
    main()
