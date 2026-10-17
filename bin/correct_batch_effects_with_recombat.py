#!/usr/bin/env python3



"""
correct_batch_effects.py: apply reCombat to count data
"""

import argparse
from pathlib import Path
import logging

from common import export_parquet
import config

import polars as pl
import numpy as np
from recombat import ReComBat

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


MODEL = 'elastic_net'
CONFIG = {
    'alpha': 1e-9,
    'copy_X': False, # to avoid unnecessary data duplication since original data are not used again
    'selection': 'random', # for faster convergence
    'random_state': 42 # for reproducible output across multiple runs
}
MAX_ITER = 1000 # RECOMBAT DEFAULT: 1000
CONV_CRITERION = 1e-4 # RECOMBAT DEFAULT: 1e-4

NB_ORPHANS_OUTFILE = "nb_excluded_orphan_samples.txt"
NB_CLIPPED_OUTFILE = "nb_clipped_samples.txt"
OUTFILE = "counts.corrected.parquet"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Parquet gathering all count data"
    )
    parser.add_argument(
        "--design", type=Path, dest="design_file", required=True, help="Design file"
    )
    parser.add_argument(
        "--out", type=Path, dest="outfile", required=True, help="Output file"
    )
    return parser.parse_args()


def remove_orphan_samples(df: pl.DataFrame, design_df: pl.DataFrame) -> tuple[pl.DataFrame, int]:
    n_samples_per_batch_df = design_df.group_by('batch').agg(n_samples=pl.col('sample').n_unique())
    batch_having_unique_samples = n_samples_per_batch_df.filter(pl.col('n_samples') == 1)['batch'].to_list()
    for batch in batch_having_unique_samples:
        batch_row = design_df.filter(pl.col('batch') == batch)
        sample = batch_row['sample'].item()
        logger.warning(f"Sample {sample} in batch {batch} is orphan. Remowing this sample.")
        df = df.drop(sample)
    return df, len(batch_having_unique_samples)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():

    args = parse_args()

    # ----------------------------------------------------
    # PARSING DATA
    # ----------------------------------------------------

    logger.info(f"Parsing counts from {args.count_file}")
    df = pl.read_parquet(args.count_file)
    logger.info(f"Getting design from {args.design_file}")
    design_df = pl.read_csv(args.design_file)

    # ----------------------------------------------------
    # REMOVING SAMPLES THAT ARE UNIQUE IN THEIR BATCH
    # ----------------------------------------------------

    df, nb_orphans = remove_orphan_samples(df, design_df)
    logger.info(f"Removed {nb_orphans} samples that were unique in their batch")

    # ----------------------------------------------------
    # ARANGING DATA
    # ----------------------------------------------------

    # filling Nan with None
    df = df.fill_nan(None)
    logger.info("Filtering out genes with all-zero counts")
    n_genes = df.height
    df = df.filter(pl.any_horizontal(pl.exclude(config.GENE_ID_COLNAME) != 0))
    logger.info(f"Removed {n_genes - df.height} genes with all-zero counts")

    # ----------------------------------------------------
    # align the design on the matrix columns (same order, same samples)
    # ----------------------------------------------------

    samples = [col for col in df.columns if col != config.GENE_ID_COLNAME]
    design_df = (
        pl.DataFrame({"sample": samples})
        .join(design_df, on="sample", how="left")
    )

    unrelated_samples_df = design_df.filter(pl.col('batch').is_null())
    if not unrelated_samples_df.is_empty():
        unrelated_samples = unrelated_samples_df["sample"].to_list()
        raise ValueError(f"Some samples were not present in the design dataframe: {unrelated_samples}.")

    # ----------------------------------------------------
    # RECOMBAT
    # ----------------------------------------------------

    # get batches in the same order as in the count dataset
    batch_series = design_df['batch']

    if len(batch_series.unique()) == 1:
        logger.info("Only one batch found. Skipping batch correction.")
        logger.info("Re-scaling data to have only positive values.")
        #df = df.with_columns([(2 ** pl.col(col) - 1).alias(col) for col in samples])

    else:

        model = ReComBat(
            model=MODEL,
            config=CONFIG,
            conv_criterion=CONV_CRITERION,
            max_iter=MAX_ITER
        )

        gene_ids = df.select(config.GENE_ID_COLNAME)
        df = df.drop(config.GENE_ID_COLNAME)

        kwargs = {
            'data': df.to_numpy().T,
            'batches': batch_series.to_numpy()
        }

        x = model.fit_transform(**kwargs)

        # ----------------------------------------------------
        # RE-SCALING DATA TO HAVE ONLY POSITIVE VALUES
        # ----------------------------------------------------

        # performing 2 ** x - 1 on all data
        logger.info("Re-scaling data to have only positive values.")
        #x = np.exp2(x) - 1

        transformed_df_without_gene_ids = pl.DataFrame(x.T, schema=df.columns)
        df = pl.concat(
            [gene_ids, transformed_df_without_gene_ids],
            how='horizontal',
            strict=True
        )

    # ----------------------------------------------------
    # CLIP NEGATIVE VALUES TO 0 (SAFETY NET)
    # ----------------------------------------------------

    n_genes = df.height
    n_samples = df.width
    n_values = n_genes * n_samples

    n_neg_values = df.select(pl.sum_horizontal(pl.col(samples) < 0)).sum().item()
    logger.info(f"Clipping {n_neg_values} ({100 * n_neg_values / n_values:.2f}%) negative values to 0")
    df = df.with_columns(pl.col(samples).clip(0))

    # ----------------------------------------------------
    # EXPORT
    # ----------------------------------------------------

    with open(NB_ORPHANS_OUTFILE, 'w') as fout:
        fout.write(str(nb_orphans))

    with open(NB_CLIPPED_OUTFILE, 'w') as fout:
        fout.write(str(n_neg_values))

    logger.info(f"Exporting corrected counts to {args.outfile}")
    export_parquet(df, args.outfile)

    logger.info("Done")


if __name__ == "__main__":
    main()
