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
from recombat import ReComBat

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#MODEL = 'linear'
MODEL = 'ridge'
#CONFIG = None
CONFIG = {'alpha': 1e-9}
MAX_ITER = 1000 # RECOMBAT DEFAULT: 1000
CONV_CRITERION = 1e-4 # RECOMBAT DEFAULT: 1e-4

OUTPUT_FILE = "counts.recombat.parquet"

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
        "--cpus", type=int, dest="nb_cpus", required=True, help="Number of CPUs"
    )
    return parser.parse_args()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():

    args = parse_args()

    df = pl.read_parquet(args.count_file)
    # filling Nan with None
    df = df.fill_nan(None)
    logger.info("Filtering out genes with all-zero counts")
    df = df.filter(pl.any_horizontal(pl.exclude(config.GENE_ID_COLNAME) != 0))

    logger.info(f"Getting design from {args.design_file}")
    design_df = pl.read_csv(args.design_file)

    # align the design on the matrix columns (same order, same samples)
    design_df = (
        pl.DataFrame({"sample": df.columns})
        .join(design_df, on="sample", how="left")
    )
    # get batches in the same order as in the count dataset
    batch_series = design_df['batch']

    if len(batch_series.unique()) == 1:
        logger.info("Only one batch found. Skipping batch correction.")
        transformed_df = df

    else:
        model = ReComBat(
            parametric=False,
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
        transformed_data = model.fit_transform(**kwargs)
        transformed_df_without_gene_ids = pl.DataFrame(transformed_data.T, schema=df.columns)
        transformed_df = pl.concat([gene_ids, transformed_df_without_gene_ids], how='horizontal', strict=True)

    export_parquet(transformed_df, OUTPUT_FILE)
    logger.info("Done")


if __name__ == "__main__":
    main()
