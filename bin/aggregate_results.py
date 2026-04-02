#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import config
import polars as pl
import yaml
from common import write_float_csv

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# outfile names
ALL_GENE_SUMMARY_OUTFILENAME = "all_genes_summary.csv"
SUMMARY_OUTFILENAME_SUFFIX = "most_stable_genes_summary.csv"
COUNTS_OUTFILENAME_SUFFIX = "most_stable_genes_transposed_counts.csv"
CUSTOM_CONTENT_MULTIQC_CONFIG_FILE = "custom_content_multiqc_config.yaml"

# quantile intervals
NB_EXPRESSION_QUANTILES = 100
NB_TOP_GENES_TO_SHOW_IN_BOX_PLOTS = 25

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
        "--target-genes",
        type=str,
        nargs="+",
        dest="target_genes",
        default=[],
        help="File containing target genes",
    )
    parser.add_argument(
        "--stats-with-scores",
        type=Path,
        nargs="+",
        dest="stat_score_files",
        required=True,
        help="Files containing statistics for all genes and stability scores by candidate genes, one per section",
    )
    parser.add_argument(
        "--multiqc-config",
        type=Path,
        dest="multiqc_config",
        required=True,
        help="MultiQC config file for custom content",
    )
    parser.add_argument(
        "--platform-stats",
        type=Path,
        dest="platform_stat_files",
        nargs="+",
        help="File containing base statistics for all genes and for all datasets for a specific platform",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        dest="metadata_files",
        help="Metadata file",
    )
    parser.add_argument(
        "--mappings", type=str, dest="mapping_files", help="Mapping file"
    )
    return parser.parse_args()


def parse_stat_score_file(file: Path) -> pl.DataFrame:
    return pl.read_csv(file).with_columns(
        pl.col(config.GENE_ID_COLNAME).cast(pl.String())
    )


def get_non_empty_dataframes(files: list[Path]) -> list[pl.DataFrame]:
    dfs = [pl.read_csv(file) for file in files]
    return [df for df in dfs if not df.is_empty()]


def cast_cols_to_string(df: pl.DataFrame) -> pl.DataFrame:
    return df.select(
        [pl.col(column).cast(pl.String) for column in df.collect_schema().names()]
    )


def concat_cast_to_string_and_drop_duplicates(files: list[Path]) -> pl.DataFrame:
    """Concatenate DataFrames, cast all columns to String, and drop duplicates.

    The first step is to concatenate the DataFrames. Then, the dataframe is cast
    to String to ensure that all columns have the same data type. Finally, duplicate
    rows are dropped.
    """
    dfs = get_non_empty_dataframes(files)
    dfs = [cast_cols_to_string(df) for df in dfs]
    concat_df = pl.concat(dfs)
    # dropping duplicates
    # casting all columns to String
    return concat_df.unique()


def cast_count_columns_to_float(df: pl.DataFrame) -> pl.DataFrame:
    return df.select(
        pl.col(config.GENE_ID_COLNAME),
        pl.exclude(config.GENE_ID_COLNAME).cast(pl.Float64),
    )


def join_data_on_gene_id(stat_df: pl.DataFrame, *dfs: pl.DataFrame) -> pl.DataFrame:
    """Merge the statistics dataframe with the metadata dataframe and the mapping dataframe."""
    # we need to ensure that the index of stat_df are strings
    for df in dfs:
        stat_df = stat_df.join(df, on=config.GENE_ID_COLNAME, how="left")
    return stat_df


def get_counts(file: Path) -> pl.DataFrame:
    # sorting dataframe (necessary to get consistent output)
    return pl.read_parquet(file).sort(config.GENE_ID_COLNAME, descending=False)


def get_metadata(metadata_files: list[Path]) -> pl.DataFrame | None:
    """Retrieve and concatenate metadata from a list of metadata files."""
    if not metadata_files:
        return None
    return concat_cast_to_string_and_drop_duplicates(metadata_files)


def get_mappings(mapping_files: list[Path]) -> pl.DataFrame | None:
    if not mapping_files:
        return None
    concat_df = concat_cast_to_string_and_drop_duplicates(mapping_files)
    # group by new gene IDs and gets the lis
    # convert the list column to a string representation
    # separate the original gene IDs with a semicolon
    return concat_df.group_by(config.GENE_ID_COLNAME).agg(
        pl.col(config.ORIGINAL_GENE_ID_COLNAME)
        .unique()
        .sort()
        .str.join(";")
        .alias(config.ORIGINAL_GENE_IDS_COLNAME)
    )


def get_status(quantile_interval: int) -> str:
    """Return the expression level status of the gene given its quantile interval."""
    if NB_EXPRESSION_QUANTILES - 5 <= quantile_interval:
        return "Very high expression"
    elif (
        NB_EXPRESSION_QUANTILES - 10 <= quantile_interval < NB_EXPRESSION_QUANTILES - 5
    ):
        return "High expression"
    elif 4 < quantile_interval <= 9:
        return "Low expression"
    elif quantile_interval <= 4:
        return "Very low expression"
    else:
        return "Medium range"


def add_expression_level_status(df: pl.DataFrame) -> pl.DataFrame:
    logger.info("Adding expression level status")
    mapping_dict = {
        quantile_interval: get_status(quantile_interval)
        for quantile_interval in range(NB_EXPRESSION_QUANTILES)
    }
    return df.with_columns(
        pl.col(config.EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
        .replace_strict(mapping_dict)
        .alias(config.EXPRESSION_LEVEL_STATUS_COLNAME)
    )


def complement_gene_summary_table(
    stat_summary_df: pl.DataFrame, *dfs: pl.DataFrame
) -> pl.DataFrame:
    """
    Add various metadata to statistics summary.
    """
    # add gene name, description and original gene IDs to statistics summary
    stat_summary_df = join_data_on_gene_id(stat_summary_df, *dfs)
    stat_summary_df = add_expression_level_status(stat_summary_df)
    return stat_summary_df


def get_most_stable_genes_counts(
    log_count_df: pl.DataFrame, stat_summary_df: pl.DataFrame
) -> pl.DataFrame:
    # getting list of top stable genes with their order
    top_genes_with_order = (
        stat_summary_df.head(NB_TOP_GENES_TO_SHOW_IN_BOX_PLOTS)
        .select(config.GENE_ID_COLNAME)
        .with_row_index("sort_order")
    )

    # join to get only existing genes and maintain order
    sorted_transposed_counts_df = log_count_df.join(
        top_genes_with_order, on=config.GENE_ID_COLNAME, how="inner"
    ).sort("sort_order", descending=False)

    # get the actual gene names that were found (in order)
    actual_gene_names = (
        sorted_transposed_counts_df.select(config.GENE_ID_COLNAME).to_series().to_list()
    )
    return sorted_transposed_counts_df.drop(
        ["sort_order", config.GENE_ID_COLNAME]
    ).transpose(column_names=actual_gene_names)


def format_multiqc_section(
    section: str, nb_sections: int, template_dict: dict, found_target_genes: list[dict]
):
    section_dict = dict(template_dict)

    parent_id = section.replace("_", " ")
    parent_name = (
        f"{section.replace('_', ' ').capitalize()} / {nb_sections}: most stable genes"
    )
    parent_description = (
        f"Most stable genes and distribution of their normalised counts for {section.replace('_', ' ')} / {nb_sections}"
        + " (section 1 corresponding to the most expressed genes)"
    )

    additional_name = ""
    if found_target_genes:
        additional_names = [
            f"{d['target_gene']} ({d['gene']})" for d in found_target_genes
        ]
        additional_name = ". Comprises " + ", ".join(additional_names)

    section_dict["parent_id"] = parent_id
    section_dict["parent_name"] = parent_name + additional_name
    section_dict["parent_description"] = parent_description

    return section_dict


def format_multiqc_sp(section: str, template_dict: dict):
    sp_dict = dict(template_dict)
    sp_dict["fn"] = sp_dict["fn"].replace("SECTION", section)
    return sp_dict


def format_genes(genes: list[str]):
    # str.maketrans("", "", "-_.") makes a mapping table for str.translate() that
    # removes all occurrences of any character in "-_." from the input string
    # it's faster than re.sub
    return pl.Series(
        [gene.lower().translate(str.maketrans("", "", "-_.")).strip() for gene in genes]
    )


def search_target_genes(df: pl.DataFrame, target_genes: list[str]) -> list[dict]:
    """
    Search for target genes in a DataFrame.

    Args:
        df (pl.DataFrame): The DataFrame to search in.
        target_genes (list[str]): The list of target genes to search for.

    Returns:
        list[dict]: A list of dictionaries associating each found target gene with its corresponding gene ID in the datasets.
    """

    unique_gene_ids = set(df[config.GENE_ID_COLNAME].to_list())

    if config.GENE_NAME_COLNAME in df.columns:
        unique_gene_ids |= set(df[config.GENE_NAME_COLNAME].to_list())

    if config.ORIGINAL_GENE_IDS_COLNAME in df.columns:
        original_gene_ids = (
            df.select(
                pl.col(config.ORIGINAL_GENE_IDS_COLNAME).str.split(by=",").explode()
            )
            .to_series()
            .to_list()
        )
        unique_gene_ids |= set(original_gene_ids)

    all_unique_gene_ids = [gene for gene in unique_gene_ids if gene is not None]

    formated_gene_ids_df = pl.DataFrame({"gene": all_unique_gene_ids}).with_columns(
        pl.col("gene")
        .map_batches(
            lambda x: format_genes(x),
            return_dtype=pl.String,
        )
        .alias("formatted_gene")
    )

    formated_target_genes_df = pl.DataFrame({"target_gene": target_genes}).with_columns(
        pl.col("target_gene")
        .map_batches(
            lambda x: format_genes(x),
            return_dtype=pl.String,
        )
        .alias("formatted_gene")
    )

    return (
        formated_gene_ids_df.join(
            formated_target_genes_df, on="formatted_gene", how="inner"
        )
        .select(["target_gene", "gene"])
        .to_dicts()
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    # --------------------------------------------------
    # Parsing counts
    # --------------------------------------------------

    count_df = get_counts(args.count_file)
    # reducing dataframe size (it is only used for plotting by MultiQC)
    count_df = cast_count_columns_to_float(count_df)

    # --------------------------------------------------
    # Parsing statistics and scores, section by section
    # --------------------------------------------------

    stat_score_dfs = []
    sections = []
    for file in args.stat_score_files:
        # the section name is at the beginning of the file name
        section = file.name.split(".")[0]
        df = parse_stat_score_file(file)
        df = df.with_columns(pl.lit(section).alias(config.SECTION_COLNAME))
        stat_score_dfs.append(df)
        sections.append(section)

    stat_score_df = pl.concat(stat_score_dfs)

    if stat_score_df.select(config.GENE_ID_COLNAME).is_duplicated().any():
        raise ValueError("Duplicate gene IDs found in statistics and scores files.")

    # sorting sections in the order (from 1 to <max nb of section>)
    sections = sorted(sections, key=lambda section: int(section.split("_")[-1]))

    # --------------------------------------------------
    # Parsing MultiQC template config for custom content
    # --------------------------------------------------

    with open(args.multiqc_config, "r") as f:
        multiqc_config = yaml.safe_load(f.read())

    # putting template parts aside
    ranking_dict = multiqc_config["custom_data"][
        "ranked_most_stable_genes_summary_template"
    ]
    ranking_sp_dict = multiqc_config["sp"]["ranked_most_stable_genes_summary_template"]
    expr_distrib_dict = multiqc_config["custom_data"][
        "expr_distrib_most_stable_genes_template"
    ]
    expr_distrib_sp_dict = multiqc_config["sp"][
        "expr_distrib_most_stable_genes_template"
    ]
    other_sections = multiqc_config["custom_content"]["order"]

    del multiqc_config["custom_data"]["ranked_most_stable_genes_summary_template"]
    del multiqc_config["sp"]["ranked_most_stable_genes_summary_template"]
    del multiqc_config["custom_data"]["expr_distrib_most_stable_genes_template"]
    del multiqc_config["sp"]["expr_distrib_most_stable_genes_template"]
    del multiqc_config["custom_content"]["order"]

    # filling dynamically the number of genes to show in box plots
    expr_distrib_dict["description"] = expr_distrib_dict["description"].replace(
        "NB_GENES", str(NB_TOP_GENES_TO_SHOW_IN_BOX_PLOTS)
    )

    # --------------------------------------------------
    # Parsing statistics per platform
    # --------------------------------------------------

    platform_datasets_stat_dfs = [
        parse_stat_score_file(file)
        for file in args.platform_stat_files
        if file is not None
    ]

    # --------------------------------------------------
    # Parsing metadata and mapping files
    # --------------------------------------------------

    metadata_files = (
        [Path(file) for file in args.metadata_files.split(" ")]
        if args.metadata_files is not None
        else []
    )
    mapping_files = (
        [Path(file) for file in args.mapping_files.split(" ")]
        if args.mapping_files is not None
        else []
    )

    # parsing metadata and mapping files
    metadata_df = get_metadata(metadata_files)
    mapping_df = get_mappings(mapping_files)
    optional_dfs = [df for df in [metadata_df, mapping_df] if df is not None]

    # --------------------------------------------------
    # Adding metadata, mapping and platform statistics information to gene summary table
    # --------------------------------------------------

    additional_data_dfs = optional_dfs + platform_datasets_stat_dfs
    all_genes_summary_df = complement_gene_summary_table(
        stat_score_df, *additional_data_dfs
    )
    print(additional_data_dfs)
    logger.info(f"Exporting statistics of all genes to: {ALL_GENE_SUMMARY_OUTFILENAME}")
    # sorting values in order to having consistent output
    all_genes_summary_df.sort(by=config.GENE_ID_COLNAME).write_csv(
        ALL_GENE_SUMMARY_OUTFILENAME, float_precision=config.CSV_FLOAT_PRECISION
    )

    # --------------------------------------------------
    # Getting summary table and counts for each section
    # Adding new sections in MultiQC config for each new expression section
    # --------------------------------------------------

    nb_sections = len(sections)
    new_mqc_config_sections = {}
    new_mqc_config_sp = {}

    logger.info("Making new sections in the MultiQC config")
    for section in sections:
        # getting best candidates for this section

        section_df = (
            all_genes_summary_df.filter(pl.col("section") == section)
            .drop("section")
            .sort(config.STABILITY_SCORE_COLNAME, nulls_last=True, maintain_order=True)
        )

        found_target_genes = []
        if args.target_genes:
            found_target_genes = search_target_genes(section_df, args.target_genes)

        section_most_stable_genes_counts_df = get_most_stable_genes_counts(
            count_df, section_df
        )

        section_summary_outfile = f"{section}.{SUMMARY_OUTFILENAME_SUFFIX}"
        write_float_csv(section_df, section_summary_outfile)

        section_counts_outfile = f"{section}.{COUNTS_OUTFILENAME_SUFFIX}"
        write_float_csv(section_most_stable_genes_counts_df, section_counts_outfile)

        # making new sections in the MultiQC config
        new_mqc_config_sections[f"genes_{section}"] = format_multiqc_section(
            section, nb_sections, ranking_dict, found_target_genes
        )
        new_mqc_config_sections[f"normalised_expr_distrib_{section}"] = (
            format_multiqc_section(
                section, nb_sections, expr_distrib_dict, found_target_genes
            )
        )
        new_mqc_config_sp[f"genes_{section}"] = format_multiqc_sp(
            section, ranking_sp_dict
        )
        new_mqc_config_sp[f"normalised_expr_distrib_{section}"] = format_multiqc_sp(
            section, expr_distrib_sp_dict
        )

    # adding new sections
    multiqc_config["custom_data"] = (
        new_mqc_config_sections | multiqc_config["custom_data"]
    )
    # specifying the filenames linked to the new sections
    multiqc_config["sp"] = new_mqc_config_sp | multiqc_config["sp"]
    # specifying the section order
    multiqc_config["custom_content"]["order"] = (
        list(new_mqc_config_sections.keys()) + other_sections
    )

    with open(CUSTOM_CONTENT_MULTIQC_CONFIG_FILE, "w") as f:
        yaml.dump(multiqc_config, f, indent=4, sort_keys=False)

    logger.info("Done")


if __name__ == "__main__":
    main()
