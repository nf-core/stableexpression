import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Get variation coefficient from count data for each gene
# The variation coefficient is the ratio of the standard deviation to the mean
# We want genes that are neither expressed too much nor too little
# Metadata (name and description) are used to annotate the results
# Likewise, mappings (original gene ids) are used to better associate gene ids with their original gene ids
# Usage:
# python get_variation_coefficients.py --counts <count_file> --metadata <metadata_file> --mapping <mapping_file>


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get variation coefficient from count data for each gene"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        dest="metadata_file",
        required=True,
        help="Metadata file",
    )
    parser.add_argument(
        "--mapping", type=Path, dest="mapping_file", required=True, help="Mapping file"
    )
    return parser.parse_args()


def get_counts(count_file: Path) -> pd.DataFrame:
    return pd.read_csv(count_file, index_col=0, header=0)


def get_metadata(metadata_file: Path) -> pd.DataFrame:
    return pd.read_csv(metadata_file, header=0)


def get_mappings(mapping_file: Path) -> pd.DataFrame:
    mapping_df = pd.read_csv(mapping_file, header=0)
    # group by new gene IDs and gets the list of distinct original gene IDs for each group
    # convert the list column to a string representation
    # separate the original gene IDs with a semicolon
    aggregated_df = (
        mapping_df.groupby("new")["original"]
        .apply(lambda x: ";".join(sorted(x.unique())))
        .reset_index()
        .rename(columns={"original": "original_gene_ids"})
    )
    return aggregated_df


def average_log2(row):
    # the dataframe has already been filtered to exclude rows where mean is 0
    return np.mean(np.log2(row + 1))  # adds 1 to avoid log(0) and to stabilize variance


def get_variation_coefficient(count_df: pd.DataFrame) -> pd.DataFrame:
    logger.info("Getting coefficients of variation")
    # handling NA values (genes that are not found in all datasets)
    # replace NaN values with 0
    count_df.fillna(0, inplace=True)

    # calculate the coefficient of variation (cv)
    # as the ratio of the standard deviation to the mean
    row_means = count_df.mean(axis=1)
    row_sds = count_df.std(axis=1)
    cv = row_sds / row_means

    # get the average log cpm value for each gene
    # to get an idea of the overall expression level of each gene
    av_log_cpm = count_df.apply(average_log2, axis=1)

    # combine results into a dataframe
    df = pd.DataFrame({"variation_coefficient": cv, "average_log_cpm": av_log_cpm})
    # order dataframe (from lowest to highest variation coefficient)
    df = df.sort_values("variation_coefficient", ascending=True)

    return df


def merge_data(
    cv_df: pd.DataFrame, metadata_df: pd.DataFrame, mapping_df: pd.DataFrame
) -> pd.DataFrame:
    # we need to ensure that the index of cv_df are strings
    return (
        cv_df.reset_index()
        .rename(columns={"index": "ensembl_gene_id"})
        .merge(metadata_df, left_on="ensembl_gene_id", right_on="gene_id", how="left")
        .merge(mapping_df, left_on="ensembl_gene_id", right_on="new", how="left")
        .drop(columns="new")
    )


def export_data(cv_df: pd.DataFrame, count_df: pd.DataFrame):
    count_outfilename = "all_normalised_counts.csv"
    logger.info(f"Exporting normalised counts to: {count_outfilename}")
    count_df.to_csv(count_outfilename, index=True)

    cv_outfilename = "variation_coefficients.csv"
    logger.info(f"Exporting variation coefficients to: {cv_outfilename}")
    cv_df.to_csv(cv_outfilename, index=False)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    count_df = get_counts(args.count_file)
    cv_df = get_variation_coefficient(count_df)

    metadata_df = get_metadata(args.metadata_file)
    mapping_df = get_mappings(args.mapping_file)

    cv_df = merge_data(cv_df, metadata_df, mapping_df)
    print(cv_df)
    export_data(cv_df, count_df)


if __name__ == "__main__":
    main()
