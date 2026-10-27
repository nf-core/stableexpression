from functools import lru_cache

import pandas as pd
import polars as pl
from src.utils import config


@lru_cache(maxsize=None)
class DataManager:
    def __init__(self):
        self.all_counts_lf: pl.LazyFrame = self.get_all_count_data()
        self.all_genes_stat_df: pl.DataFrame = self.get_all_genes_stat_data()

    @staticmethod
    def get_all_count_data() -> pl.LazyFrame:
        file = f"{config.DATA_FOLDER}/{config.ALL_COUNT_FILENAME}"
        return pl.scan_parquet(file)

    def get_sorted_samples(self) -> list[str]:
        return sorted(
            self.all_counts_lf.select(pl.exclude(config.GENE_ID_COLNAME))
            .collect_schema()
            .names()
        )

    def get_all_genes_stat_data(self) -> pl.DataFrame:
        file = f"{config.DATA_FOLDER}/{config.ALL_GENES_STAT_FILENAME}"
        stat_df = pl.read_csv(file)
        cols_to_select = ["rank"] + [
            col for col in stat_df.columns if col not in ["rank", "is_candidate"]
        ]
        return stat_df.select(cols_to_select)

    """
    def get_samples_grouped_by_dataset(self) -> list[dict]:

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
    """

    def get_sorted_genes(self) -> list[str]:
        return (
            self.all_genes_stat_df.sort(
                by=[config.RANK_COLNAME, config.SECTION_COLNAME],
                descending=False,
            )
            .select(config.GENE_ID_COLNAME)
            .to_series()
            .to_list()
        )

    def get_gene_counts(self, gene: str) -> pd.Series:
        return (
            self.all_counts_lf.filter(pl.col(config.GENE_ID_COLNAME) == gene)
            .select(pl.exclude(config.GENE_ID_COLNAME))
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

    def get_nb_sections(self) -> int:
        return self.all_genes_stat_df.select(config.SECTION_COLNAME).n_unique()

    def get_table_raw_data(self) -> list[dict]:
        return self.all_genes_stat_df.sort(
            by=[config.RANK_COLNAME, config.SECTION_COLNAME],
            descending=False,
        ).to_dicts()
