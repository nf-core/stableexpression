#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
"""
genorm.py - Compute M-measure for each gene
"""

import argparse
import logging
from pathlib import Path

import config
from common import sink_csv_with_floats, get_count_columns
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

STD_DDOF = 0
JOIN_SUFFIX = "_other"


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PARSE ARGS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--file",
        type=Path,
        dest="count_file",
        required=True,
        help="File where each row contains counts for two genes",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="outfile",
        required=True,
        help="M-measure output file",
    )
    return parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAKE PAIRS (CROSS JOIN WITHOUT DUPLICATES)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_pairs(lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Make new dataframe where each row contains the counts (across samples) of each unique pair of genes.
    The output is like:
        ┌───────────┬───────────────┬──────────────────┐
        │ gene_id   ┆ sample_A      ┆ sample_A_other   │
        │ str       ┆ ---           ┆ ---              │
        │           ┆ f32           ┆ f32              │
        ╞═══════════╪═══════════════╪══════════════════│
        │ AT1G01700 ┆ 2.057693      ┆ 2.213529         │
        ...
    """
    lf = lf.with_row_index("_idx")
    return (
        lf.join_where(
            lf,
            pl.col("_idx") < pl.col("_idx_other"),
            suffix=JOIN_SUFFIX
        )
        .drop("_idx", "_idx_other")
    )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMPUTATION OF LOGS OF PAIRWISE EXPRESSION RATIOS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def compute_ratios(cross_join_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Perform a cross join of a count file with itself and compute log ratio for each pair of gene, for each column

    The output is like:
        gene_id	gene_id_other	sampleA	sampleB
        AT1G01060	AT4G34290	0.48072526	-0.049063634
        AT1G01060	AT4G34460	0.21205564	0.15575644
    """
    column_pairs = {
        col: f"{col}{JOIN_SUFFIX}"
        for col in get_count_columns(cross_join_lf)
        if not col.endswith(JOIN_SUFFIX)
    }
    return (
        cross_join_lf.with_columns( # computing log for each count column
            pl.exclude([config.GENE_ID_COLNAME, f"{config.GENE_ID_COLNAME}{JOIN_SUFFIX}"]).log(base=2)
        ).select(
            [pl.col(config.GENE_ID_COLNAME), pl.col(f"{config.GENE_ID_COLNAME}{JOIN_SUFFIX}")]
            + [
                (pl.col(col) - pl.col(other_col)).alias(f"{col}_log_ratio")
                for col, other_col in column_pairs.items()
            ]
        )
    )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMPUTATION OF STANDARD DEVIATION OF RATIOS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def compute_standard_deviations(ratios_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    For each gene pair, combine all the *_log_ratio columns into one list and computes the standard deviation across samples.
    The result is one row per pair: (ratios_std, gene_id, gene_id_other)
    Computing the std manually for better performance with polars (as long as pl.std_horizontal does not exist)

    The output is like:
        gene_id gene_id_other ratio_stds
        AT1G01060   AT1G01061   [0.04122424, 0.049532797, 0.050319236, ...]
    """
    ratio_columns = [
        c for c in ratios_lf.collect_schema().names() if c.endswith("_log_ratio")
    ]
    if not ratio_columns:
        raise ValueError("No '*_log_ratio' columns found in input")

    n = len(ratio_columns)

    mean = pl.mean_horizontal(ratio_columns)
    variance = pl.sum_horizontal([(pl.col(c) - mean) ** 2 for c in ratio_columns]) / (n - STD_DDOF)

    return ratios_lf.select(
        pl.col(config.GENE_ID_COLNAME),
        pl.col(f"{config.GENE_ID_COLNAME}_other"),
        variance.sqrt().cast(pl.Float64).alias(config.RATIOS_STD_COLNAME),
    )




def compute_m_measures(std_lf: pl.LazyFrame, n_genes: int) -> pl.LazyFrame:
    """
    Compute the m-measure for each gene.
    The m-measure is the sum of the ratios standard deviations divided by the number of ratios minus 1.
    The output is like:
        ┌───────────┬──────────────────┐
        │ gene_id   ┆ genorm_m_measure │
        │ ---       ┆ ---              │
        │ str       ┆ f32              │
        ╞═══════════╪══════════════════╡
        │ AT1G01380 ┆ 0.000144         │
        │ AT1G01700 ┆ 0.000154         │
        │ AT1G01750 ┆ 0.000154         │
        ...
    """

    std = pl.col(config.RATIOS_STD_COLNAME)

    both_sides = pl.concat(
        [
            std_lf.select(pl.col(config.GENE_ID_COLNAME), std),
            std_lf.select(pl.col(f"{config.GENE_ID_COLNAME}{JOIN_SUFFIX}").alias(config.GENE_ID_COLNAME), std),
        ]
    )

    return (
        both_sides.group_by(config.GENE_ID_COLNAME)
        .agg((std.sum() / (n_genes - 1)).alias(config.GENORM_M_MEASURE_COLNAME))
        .sort(config.GENE_ID_COLNAME)
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info(f"Computing M-measure for {args.count_file}")

    count_lf = pl.scan_parquet(args.count_file)
    n_genes = count_lf.select(pl.len()).collect().item()

    # replace 0 values by a very small count value
    count_lf = count_lf.with_columns(pl.exclude("gene_id").replace(0, 1e-8))

    unique_pairs_lf = make_pairs(count_lf)

    ratios_lf = compute_ratios(unique_pairs_lf)

    ratio_std_lf = compute_standard_deviations(ratios_lf)

    # computing M measures for these gene IDs
    m_measure_lf = compute_m_measures(ratio_std_lf, n_genes)
    print(m_measure_lf.collect())

    sink_csv_with_floats(m_measure_lf, args.outfile, float_precision=5)

    logger.info("Done")


if __name__ == "__main__":
    main()
