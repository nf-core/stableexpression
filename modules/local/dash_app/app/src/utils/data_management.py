import polars as pl
import pandas as pd
from functools import lru_cache

from src.utils import config


@lru_cache(maxsize=None)
class DataManager:
    def __init__(self):
        self.all_counts_lf: pl.LazyFrame = self.get_all_count_data()
        self.grouped_samples: list[dict] = self.get_samples_grouped_by_dataset()
        self.candidate_genes_stat_df: pl.DataFrame = (
            self.get_candidate_genes_stat_data()
        )
        self.all_gene_stats_df: pl.DataFrame = self.get_all_genes_stat_data()
        self.genes: list[str] = self.get_sorted_genes()

    @staticmethod
    def get_all_count_data() -> pl.LazyFrame:
        file = f"{config.DATA_FOLDER}/{config.ALL_COUNT_FILENAME}"
        return pl.scan_parquet(file)

    def get_samples_in_count_data(self) -> list[str]:
        return (
            self.all_counts_lf.select(pl.exclude(config.ENSEMBL_GENE_ID_COLNAME))
            .collect_schema()
            .names()
        )

    def get_candidate_genes_stat_data(self) -> pl.DataFrame:
        file = f"{config.DATA_FOLDER}/{config.CANDIDATE_GENES_STAT_FILENAME}"
        stat_df = pl.read_csv(file)
        cols_to_select = ["rank"] + [
            col for col in stat_df.columns if col not in ["rank", "is_candidate"]
        ]
        return stat_df.select(cols_to_select)

    def get_all_genes_stat_data(self) -> pl.DataFrame:
        file = f"{config.DATA_FOLDER}/{config.ALL_GENES_STAT_FILENAME}"
        return pl.read_csv(file)

    def get_samples_grouped_by_dataset(self) -> list[dict]:
        samples_in_count_data = self.get_samples_in_count_data()
        samples_grouped_by_dataset = []

        design_file = f"{config.DATA_FOLDER}/{config.ALL_DESIGNS_FILENAME}"
        design_df = pd.read_csv(design_file)

        for group, samples in design_df.groupby(["batch", "condition"])["sample"]:
            batch, condition = group  # unpacking
            batch_condition_samples_dict = {
                "group": f"Dataset: {batch} || Condition: {condition}",
                "items": [
                    {"value": sample, "label": sample}
                    for sample in samples.to_list()
                    if sample in samples_in_count_data
                ],
            }
            samples_grouped_by_dataset.append(batch_condition_samples_dict)

        return samples_grouped_by_dataset

    def get_sorted_genes(self) -> list[str]:
        return (
            self.candidate_genes_stat_df.sort(
                by=config.STABILITY_SCORE_COLNAME, descending=False
            )
            .select(config.ENSEMBL_GENE_ID_COLNAME)
            .to_series()
            .to_list()
        )

    def get_gene_counts(self, gene: str) -> pd.Series:
        return (
            self.all_counts_lf.filter(pl.col(config.ENSEMBL_GENE_ID_COLNAME) == gene)
            .collect()
            .to_pandas()
            .iloc[0]
        )

    def get_sample_counts(self, sample: str) -> pd.Series:
        return (
            self.all_counts_lf.select(sample)
            .drop_nulls()
            .collect()
            .to_pandas()
            .iloc[:, 0]
        )
