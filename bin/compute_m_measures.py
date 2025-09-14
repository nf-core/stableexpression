#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import polars as pl
from pathlib import Path
import argparse
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
RATIOS_STD_COLNAME = "ratios_stds"
M_MEASURE_COLNAME = "m_measure"

M_MEASURE_OUTFILE_NAME = "m_measures.csv"

DEFAULT_CHUNKSIZE = 300
NB_GENE_ID_CHUNK_FOLDERS = 100


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
        dest="count_file",
        required=True,
        help="File containing std of lof expression ratios",
    )
    parser.add_argument(
        "--std-files",
        type=str,
        dest="std_files",
        required=True,
        help="File containing std of lof expression ratios",
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


def concat_all_std_data(files: list[Path], low_memory: bool) -> pl.LazyFrame:
    lfs = [pl.scan_parquet(file, low_memory=low_memory) for file in files]
    lf = pl.concat(lfs)
    return (
        lf.explode(RATIOS_STD_COLNAME)
        .group_by(ENSEMBL_GENE_ID_COLNAME)
        .agg(pl.col(RATIOS_STD_COLNAME))
    )


def compute_m_measures(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        (
            pl.col(RATIOS_STD_COLNAME).list.sum()
            / (pl.col(RATIOS_STD_COLNAME).list.len() - 1)
        ).alias(M_MEASURE_COLNAME),
    )


def get_chunks(lst: list, chunksize: int):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), chunksize):
        yield lst[i : i + chunksize]


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    low_memory = True if args.task_attempts > 1 else False
    files = [Path(file) for file in args.std_files.split(" ")]

    logger.info("Getting list of gene IDs")
    count_lf = pl.scan_parquet(args.count_file, low_memory=low_memory)

    #############################################################################
    # MAKING A FOLDER FOR EACH CHUNK OF GENE IDS
    #############################################################################
    gene_ids = count_lf.select(ENSEMBL_GENE_ID_COLNAME).collect().to_series().to_list()
    gene_ids = sorted(gene_ids)

    chunksize = max(
        1, int(len(gene_ids) / NB_GENE_ID_CHUNK_FOLDERS)
    )  # 1 if len(gene_ids) < NB_GENE_ID_CHUNK_FOLDERS
    gene_id_list_chunks = list(get_chunks(gene_ids, chunksize=chunksize))

    gene_id_chunk_folders = []
    for i in range(len(gene_id_list_chunks)):
        gene_id_chunk_folder = Path(f"gene_ids_{i}")
        gene_id_chunk_folder.mkdir(exist_ok=True)
        gene_id_chunk_folders.append(gene_id_chunk_folder)

    #############################################################################
    # EXPORTING GENE DATA TO THEIR RESPECTIVE CHUNK FOLDER
    #############################################################################
    # progressively decreasing the chunksize if OOM
    chunksize = int(DEFAULT_CHUNKSIZE / args.task_attempts)
    chunk_files_list = [
        files[i : i + chunksize] for i in range(0, len(files), chunksize)
    ]

    logger.info("Parsing std data by chunks")
    for i, chunk_files in enumerate(chunk_files_list):
        # parsing files and making a first list concatenation
        concat_lf = concat_all_std_data(chunk_files, low_memory)

        # looping through each group of gene IDs
        for j, (gene_id_list_chunk, gene_id_chunk_folder) in enumerate(
            zip(gene_id_list_chunks, gene_id_chunk_folders)
        ):
            # writing all data corresponding to this group of gene IDs in a specific folder
            outfile = gene_id_chunk_folder / f"chunk.{i}.parquet"
            concat_df = concat_lf.filter(
                pl.col(ENSEMBL_GENE_ID_COLNAME).is_in(gene_id_list_chunk)
            ).collect()
            concat_df.write_parquet(outfile)

    #############################################################################
    # GATHERING ALL DATA CHUNK BY CHUNK AND COMPUTING M MEASURE FOR EACH GENE
    #############################################################################
    computed_genes = 0
    nb_ratios_per_gene = set()
    logger.info(
        "Concatenating all std data by chunk of gene IDs and computing M measures"
    )
    with open(M_MEASURE_OUTFILE_NAME, "a") as fout:
        for i, gene_id_chunk_folder in enumerate(gene_id_chunk_folders):
            chunk_files = list(gene_id_chunk_folder.iterdir())

            concat_lf = concat_all_std_data(chunk_files, low_memory).sort(
                ENSEMBL_GENE_ID_COLNAME
            )

            # computing M measures for these gene IDs
            m_measure_lf = compute_m_measures(concat_lf)
            m_measure_df = m_measure_lf.collect()

            #################################################
            # checks
            #################################################
            if m_measure_df[ENSEMBL_GENE_ID_COLNAME].is_duplicated().any():
                raise ValueError("Duplicate values found for gene IDs!")

            process_gene_ids = sorted(
                m_measure_df.select(ENSEMBL_GENE_ID_COLNAME).to_series().to_list()
            )
            if process_gene_ids != gene_id_list_chunks[i]:
                raise ValueError("Incorrect gene IDs found!")

            computed_genes += len(m_measure_df)

            unique_nb_ratios = (
                concat_lf.with_columns(
                    pl.col(RATIOS_STD_COLNAME).list.len().alias("length")
                )
                .select("length")
                .unique()
                .collect()
                .to_series()
                .to_list()
            )
            nb_ratios_per_gene.update(unique_nb_ratios)

            #################################################
            #################################################

            # appending to output file
            if i == 0:
                m_measure_df.write_csv(fout, include_header=True)
            else:
                m_measure_df.write_csv(fout, include_header=False)

    logger.info(f"Number of gene IDs: {len(gene_ids)}")
    logger.info(f"Number of computed genes: {computed_genes}")
    if computed_genes != len(gene_ids):
        raise ValueError(
            f"Number of computed genes: {computed_genes} != number of gene IDs: {len(gene_ids)}"
        )

    if len(nb_ratios_per_gene) > 1:
        logger.warning(
            f"Got multiple number of std ratios to compute: {list(nb_ratios_per_gene)}"
        )


if __name__ == "__main__":
    main()
