PLOTLY_APP_PORT = 8080
HOST = "0.0.0.0"

LOGO_FILENAME = "assets/nf-core-stableexpression_logo_light_small.png"

LOGGING_FORMAT = "[%(asctime)s] [%(name)s] %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d_%H-%M-%S"

APP_TITLE = "Counts"
UPDATE_TITLE = "Updating ..."

DATA_FOLDER = "data"

ALL_COUNT_FILENAME = "all_counts.parquet"
CANDIDATE_GENES_STAT_FILENAME = "stats_with_scores.csv"
ALL_GENES_STAT_FILENAME = "stats_all_genes.csv"
ALL_DESIGNS_FILENAME = "whole_design.csv"

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
STD_COLNAME = "standard_deviation"
STABILITY_SCORE_COLNAME = "stability_score"

AG_GRID_DEFAULT_COLUMN_DEF = {
    "filter": True,
    "resizable": True,
    "editable": False,
    "sortable": True,
}

AG_GRID_DEFAULT_OPTIONS = {"pagination": True, "paginationAutoPageSize": True}
