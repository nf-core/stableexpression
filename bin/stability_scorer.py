#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

from typing import ClassVar
import polars as pl
from dataclasses import dataclass, field
import logging

logger = logging.getLogger(__name__)

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
GENE_COUNT_COLNAME = "count"
SAMPLE_COLNAME = "sample"


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME)).collect_schema().names()



@dataclass
class StabilityScorer:

    STAT_COLS: ClassVar[dict] = dict(
        VAR_COEFF="variation_coefficient",
        STD="standard_deviation",
        MEAN="mean",
        EXPRESSION_LEVEL_QUANTILE_INTERVAL="expression_level_quantile_interval",
        EXPRESSION_LEVEL_STATUS="expression_level_status",
        NB_NULLS="total_nb_nulls",
        NB_NULLS_VALID_SAMPLES="nb_nulls_valid_samples",
        NB_ZEROS="nb_zeros",
        STABILITY_SCORE="stability_score"
    )

    # we want to select samples that show a particularly low nb of genes
    MIN_RATIO_GENE_COUNT_TO_MEAN: ClassVar[float] = 0.75  # experimentally chosen
    WEIGHT_RATIO_NB_NULLS: ClassVar[float] = 1
    # quantile intervals
    NB_QUANTILES: ClassVar[int] = 100

    count_lf: pl.LazyFrame
    platform: str | None = field(default=None)

    gene_count_per_sample_df: pl.DataFrame = field(init=False)
    stat_lf: pl.LazyFrame = field(init=False)
    count_columns: list[str] = field(init=False)
    samples_with_low_gene_count: list[str] = field(init=False)

    def __post_init__(self):
        self.count_columns = get_count_columns(self.count_lf)
        self.gene_count_per_sample_df = self.get_gene_counts_per_sample()
        self.samples_with_low_gene_count = self.get_samples_with_low_gene_count()


    def get_colname(self, key: str) -> str:
        return f"{self.platform}_{self.STAT_COLS[key]}" if self.platform else self.STAT_COLS[key]

    def get_valid_counts(self) -> pl.LazyFrame:
        return self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME))

    def get_gene_counts_per_sample(self) -> pl.DataFrame:
        """
        Get the number of non-null values per sample.
        :return:
        A polars dataframe containing 2 columns:
            - sample: name of the sample
            - nb_not_nulls: number of non-null values
        """
        return (
            self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME))
            .count()
            .collect()
            .transpose(
                include_header=True,
                header_name="sample",
                column_names=["count"]
            )
        )

    def get_samples_with_low_gene_count(self) -> list[str]:
        mean_gene_count = self.gene_count_per_sample_df[GENE_COUNT_COLNAME].mean()
        return (
            self.gene_count_per_sample_df.filter(
                (pl.col(GENE_COUNT_COLNAME) / mean_gene_count)
                < self.MIN_RATIO_GENE_COUNT_TO_MEAN
            )
            .select(SAMPLE_COLNAME)
            .to_series()
            .to_list()
        )

    def get_main_statistics(self) -> pl.LazyFrame:
        """
        Compute count descriptive statistics for each gene in the count dataframe.
        """
        logger.info("Getting descriptive statistics")
        # computing main stats
        augmented_count_lf = self.count_lf.with_columns(
            mean=pl.concat_list(self.count_columns).list.drop_nulls().list.mean(),
            std=pl.concat_list(self.count_columns).list.drop_nulls().list.std(),
        )
        return augmented_count_lf.select(
            pl.col(ENSEMBL_GENE_ID_COLNAME),
            pl.col("mean").alias(self.get_colname("MEAN")),
            pl.col("std").alias(self.get_colname("STD")),
            (pl.col("std") / pl.col("mean")).alias(self.get_colname("VAR_COEFF")),
        )

    def compute_nb_null_values(self):
        # the samples showing a low gene count will not be taken into account for the zero count penalty
        cols_to_exclude = [ENSEMBL_GENE_ID_COLNAME] + self.samples_with_low_gene_count
        total_nb_nulls = (
            self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME).is_null())
            .collect()
            .sum_horizontal()
        )
        nb_nulls_valid_samples = (
            self.count_lf.select(pl.exclude(cols_to_exclude).is_null())
            .collect()
            .sum_horizontal()
        )
        self.stat_lf = self.stat_lf.with_columns(
            total_nb_nulls.alias(self.get_colname("NB_NULLS")),
            nb_nulls_valid_samples.alias(self.get_colname("NB_NULLS_VALID_SAMPLES")),
        )

    def get_quantile_intervals(self):
        """
        Compute the quantile intervals for the mean expression levels of each gene in the dataframe.

        The function assigns to each gene a quantile interval of its mean cpm compared to all genes.
        """
        logger.info("Getting cpm quantiles")
        self.stat_lf = self.stat_lf.with_columns(
            (pl.col(self.get_colname("MEAN")).rank() / pl.col(self.get_colname("MEAN")).count() * self.NB_QUANTILES)
            .floor()
            .cast(pl.Int8)
            # we want the only value = NB_QUANTILES to be NB_QUANTILES - 1
            # because the last quantile interval is [NB_QUANTILES - 1, NB_QUANTILES]
            .replace({self.NB_QUANTILES: self.NB_QUANTILES - 1})
            .alias(self.get_colname("EXPRESSION_LEVEL_QUANTILE_INTERVAL"))
        )

    def compute_stability_score(self):
        logger.info("Computing stability score")
        nb_valid_samples = self.gene_count_per_sample_df.select(pl.len()).item() - len(
            self.samples_with_low_gene_count
        )
        ratio_nb_nulls = (
            self.stat_lf.select(
                pl.col(self.get_colname("NB_NULLS_VALID_SAMPLES")) / nb_valid_samples
            )
            .collect()
            .to_series()
        )
        expr = (
            pl.col(self.get_colname("STD")) + ratio_nb_nulls * self.WEIGHT_RATIO_NB_NULLS
        )
        self.stat_lf = self.stat_lf.with_columns(expr.alias(self.get_colname("STABILITY_SCORE")))

    def compute_statistics_and_score(self) -> pl.LazyFrame:
        logger.info("Computing statistics and stability score")
        # getting expression statistics
        self.stat_lf = self.get_main_statistics()
        # adding column for nb of null values for each gene
        self.compute_nb_null_values()
        # computing stability score
        self.compute_stability_score()
        # getting quantile intervals
        self.get_quantile_intervals()

        return self.stat_lf
