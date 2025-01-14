#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import polars as pl
from tqdm import tqdm
import psutil
from pathlib import Path
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
M_MEASURE_OUTFILE_NAME = "m_measures.csv"

# experimentally chosen
FREE_MEMORY_TO_GENE_CHUNK_SIZE_RATIO = 400000
FREE_MEMORY_TO_EXPERSSION_RATIO_CHUNK_SIZE_RATIO = 4.5


# prepare folders
count_summary_chunk_folder = Path().cwd() / "chunks"
cross_join_chunk_folder = Path().cwd() / "cross_joins"
ratios_chunk_folder = Path().cwd() / "ratios"
std_chunk_folder = Path().cwd() / "std"

count_summary_chunk_folder.mkdir(exist_ok=True)
cross_join_chunk_folder.mkdir(exist_ok=True)
ratios_chunk_folder.mkdir(exist_ok=True)
std_chunk_folder.mkdir(exist_ok=True)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Compute M-measure for each gene")
    parser.add_argument(
        "--counts",
        type=Path,
        dest="count_summary_file",
        required=True,
        help="File containing normalised counts for all genes and all samples",
    )
    parser.add_argument(
        "--task-attempts",
        dest="task_attempts",
        type=int,
        default=1,
        help="Number of task attempts",
    )
    return parser.parse_args()


def get_nb_rows(lf: pl.LazyFrame):
    return lf.select(pl.len()).collect().item()


def get_nb_cols(lf: pl.LazyFrame):
    return len(lf.collect_schema().names())


def get_gene_chunk_size(lf: pl.LazyFrame, gene_chunk_size_factor: float):
    free_memory = psutil.virtual_memory().free
    # getting gene chunk size as a function of free memory, nb of rows and nb of cols
    ratio = free_memory / ((get_nb_rows(lf) ** 2) * (2 * get_nb_cols(lf)))
    return int(FREE_MEMORY_TO_GENE_CHUNK_SIZE_RATIO * gene_chunk_size_factor * ratio)


def get_ratio_chunk_size(lf: pl.LazyFrame, ratio_chunk_size_factor: float):
    free_memory = psutil.virtual_memory().free
    # getting gene chunk size as a function of free memory, nb of rows and nb of cols
    ratio = free_memory / (get_nb_rows(lf) * (2 * get_nb_cols(lf)))
    return int(
        FREE_MEMORY_TO_EXPERSSION_RATIO_CHUNK_SIZE_RATIO
        * ratio_chunk_size_factor
        * ratio
    )


def parse_count_summary_dataset(file: Path, low_memory: bool) -> pl.LazyFrame:
    lf = pl.scan_parquet(file, low_memory=low_memory).fill_null(0).fill_nan(0)
    count_columns = get_count_columns(lf)
    cols = [pl.col(ENSEMBL_GENE_ID_COLNAME)] + [
        pl.col(column).replace({0: 1e-8}).cast(pl.Float32) for column in count_columns
    ]
    return lf.select(cols)


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return [
        col
        for col in lf.collect_schema().names()
        if not col.startswith(ENSEMBL_GENE_ID_COLNAME)
    ]


def parse_count_summary_chunk_dataset(file: Path, low_memory: bool) -> pl.LazyFrame:
    return pl.scan_parquet(file, low_memory=low_memory)


def make_cross_join(lf1: pl.LazyFrame, lf2: pl.LazyFrame) -> pl.LazyFrame:
    return lf1.join(
        lf2, how="cross", suffix="_other"
    )  # Perform a cross join with itself


def compute_ratios(file: Path, low_memory: bool) -> pl.LazyFrame:
    # getting ratios for each sample
    cross_join_lf = pl.scan_parquet(file, low_memory=low_memory)
    column_pairs = {
        col: f"{col}_other"
        for col in get_count_columns(cross_join_lf)
        if not col.endswith("_other")
    }
    return cross_join_lf.select(
        [pl.col(ENSEMBL_GENE_ID_COLNAME), pl.col(f"{ENSEMBL_GENE_ID_COLNAME}_other")]
        + [
            (pl.col(col) / pl.col(other_col)).log(base=2).alias(f"{col}_log_ratio")
            for col, other_col in column_pairs.items()
        ]
    )


def check_has_nan_values(df: pl.DataFrame, col: str):
    df = df.fill_nan(None)
    if df.select(pl.col(col)).null_count().item() > 0:
        raise ValueError


def compute_standard_deviations(
    file: Path, low_memory: bool, ratio_chunk_size: int
) -> pl.LazyFrame:
    ratios_lf = pl.scan_parquet(file, low_memory=low_memory)
    ratio_columns = [
        col for col in ratios_lf.collect_schema().names() if col.endswith("_log_ratio")
    ]
    concat_ratios_lf = ratios_lf.select(
        [
            pl.concat_list(
                [pl.col(col) for col in ratio_columns[i : i + ratio_chunk_size]]
            ).alias(f"concat_list_chunk_{i // ratio_chunk_size}")
            for i in range(0, len(ratio_columns), ratio_chunk_size)
        ]
    ).select(pl.concat_list(pl.all()).alias("ratios"))
    return pl.concat(
        [
            concat_ratios_lf.select("ratios"),
            ratios_lf.select(
                pl.exclude("^.*_log_ratio$")
            ),  # ensembl_gene_id & ensembl_gene_id_other
        ],
        how="horizontal",
    ).select(
        pl.col("ratios").list.std().alias("ratios_std"),
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        pl.col(f"{ENSEMBL_GENE_ID_COLNAME}_other"),
    )


def group_standard_deviations(std_lf: pl.LazyFrame) -> pl.LazyFrame:
    # getting the standard devs for "a" genes
    std_a = (
        std_lf.group_by(ENSEMBL_GENE_ID_COLNAME)
        .agg("ratios_std")  # getting list of ratio std for this gene
        .select(pl.col(ENSEMBL_GENE_ID_COLNAME), pl.col("ratios_std"))
    )
    # getting the standard devs for "b" genes
    std_b = (
        std_lf.group_by(f"{ENSEMBL_GENE_ID_COLNAME}_other")
        .agg("ratios_std")  # getting list of ratio std for this gene
        .select(
            pl.col(f"{ENSEMBL_GENE_ID_COLNAME}_other").alias(ENSEMBL_GENE_ID_COLNAME),
            pl.col("ratios_std"),
        )
    )
    # concatenating both dataframes vertically
    return pl.concat([std_a, std_b], how="vertical")


def concat_all_std_data() -> pl.LazyFrame:
    lfs = [pl.scan_parquet(file) for file in list(std_chunk_folder.iterdir())]
    lf = pl.concat(lfs)
    return (
        lf.explode("ratios_std")
        .group_by(ENSEMBL_GENE_ID_COLNAME)
        .agg(pl.col("ratios_std"))
    )


def compute_m_measures(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        (
            pl.col("ratios_std").list.mean() / (pl.col("ratios_std").list.len() - 1)
        ).alias("m_measure"),
    )


#####################################################
#####################################################
# STEP FUNCTIONS
#####################################################
#####################################################


def split_count_summary_in_chunks(lf: pl.LazyFrame, gene_chunk_size: int):
    lf = lf.with_row_index(name="index")
    nb_rows = get_nb_rows(lf)
    logger.info(f"Number of rows (genes) in count file: {nb_rows}")
    nb_chunks = nb_rows // gene_chunk_size + 1
    logger.info(f"Number of chunks: {nb_chunks}")

    for i, start in tqdm(
        enumerate(range(0, nb_rows + 1, gene_chunk_size)), total=nb_chunks
    ):
        partition = (
            lf.filter(
                (pl.col("index") >= start) & (pl.col("index") < start + gene_chunk_size)
            )
            .drop("index")
            .collect()
        )
        outfile = count_summary_chunk_folder / f"count_summary_chunk_{i}.parquet"
        partition.write_parquet(outfile)


def compute_and_export_cross_join_data(low_memory: bool):
    count_summary_chunk_files = list(count_summary_chunk_folder.iterdir())
    nb_files = len(count_summary_chunk_files)
    total = int(nb_files * (nb_files - 1) / 2)

    with tqdm(total=total) as pbar:
        for i, file in enumerate(count_summary_chunk_files):
            for j, file_other in enumerate(count_summary_chunk_files):
                # Filter to keep only unique pairs (i < j)
                if i > j:
                    continue

                # parsing
                lf = parse_count_summary_chunk_dataset(file, low_memory)
                lf_other = parse_count_summary_chunk_dataset(file_other, low_memory)

                # cross join of the two dataframes
                cross_join_lf = make_cross_join(lf, lf_other)

                cross_join_df = cross_join_lf.collect()
                if len(cross_join_df) == 0:
                    raise ValueError(
                        f"No output following treatment of file {str(file)}"
                    )

                outfile = cross_join_chunk_folder / f"cross_join_{i}_{j}.parquet"
                cross_join_df.write_parquet(outfile)

                # adding 1 to the progress bar
                pbar.update()


def compute_and_export_count_ratios(low_memory: bool):
    cross_join_chunk_files = list(cross_join_chunk_folder.iterdir())
    nb_files = len(cross_join_chunk_files)

    for file in tqdm(cross_join_chunk_files, total=nb_files):
        ratios_lf = compute_ratios(file, low_memory)

        ratios_df = ratios_lf.collect()
        if len(ratios_df) == 0:
            raise ValueError(f"No output following treatment of file {str(file)}")

        outfile = ratios_chunk_folder / file.name.replace("cross_join", "ratios")
        ratios_df.write_parquet(outfile)


def compute_and_export_standard_deviations(low_memory: bool, ratio_chunk_size: int):
    ratios_chunk_files = list(ratios_chunk_folder.iterdir())
    nb_files = len(ratios_chunk_files)

    for file in tqdm(ratios_chunk_files, total=nb_files):
        std_lf = compute_standard_deviations(file, low_memory, ratio_chunk_size)
        std_lf = group_standard_deviations(std_lf)

        std_df = std_lf.collect()
        outfile = std_chunk_folder / file.name.replace("ratios", "std")
        std_df.write_parquet(outfile)


def compute_and_export_m_measure(low_memory: bool):
    all_std_lf = concat_all_std_data()
    m_measure_lf = compute_m_measures(all_std_lf)
    m_measure_lf.collect().write_csv(M_MEASURE_OUTFILE_NAME)
    print(m_measure_lf.collect())


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    low_memory = True if args.task_attempts > 1 else False

    logger.info("Parsing count file")
    lf = parse_count_summary_dataset(args.count_summary_file, low_memory)

    gene_chunk_size_factor = args.task_attempts
    ratio_chunk_size_factor = args.task_attempts

    gene_chunk_size = get_gene_chunk_size(lf, gene_chunk_size_factor)
    ratio_chunk_size = get_ratio_chunk_size(lf, ratio_chunk_size_factor)
    logger.info(f"Gene chunk size: {gene_chunk_size}")
    logger.info(f"Ratio chunk size: {ratio_chunk_size}")

    logger.info("Splitting count file into chunks")
    split_count_summary_in_chunks(lf, gene_chunk_size)

    logger.info("Computing cross join data")
    compute_and_export_cross_join_data(low_memory)

    logger.info("Computing pairwise expression ratios")
    compute_and_export_count_ratios(low_memory)

    logger.info("Computing standard deviation of pairwise expression ratios")
    compute_and_export_standard_deviations(low_memory, ratio_chunk_size)

    logger.info("Computing m-measure for all genes")
    compute_and_export_m_measure(low_memory)
    logger.info("M-measure computation done")


if __name__ == "__main__":
    main()
