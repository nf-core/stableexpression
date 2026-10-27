PLOTLY_APP_PORT = 8080
HOST = "0.0.0.0"

LOGO_FILENAME = "assets/nf-core-stableexpression_logo_light_small.png"

LOGGING_FORMAT = "[%(asctime)s] [%(name)s] %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d_%H-%M-%S"

APP_TITLE = "Counts"
UPDATE_TITLE = "Updating ..."

DATA_FOLDER = "data"

ALL_COUNT_FILENAME = "all_counts.imputed.parquet"
ALL_GENES_STAT_FILENAME = "all_genes_summary.csv"
ALL_DESIGNS_FILENAME = "whole_design.csv"

GENE_ID_COLNAME = "gene_id"
STD_COLNAME = "standard_deviation"
STABILITY_SCORE_COLNAME = "stability_score"
RANK_COLNAME = "rank"
SECTION_COLNAME = "section"

AG_GRID_DEFAULT_COLUMN_DEF = {
    "filter": True,
    "resizable": True,
    "editable": False,
    "sortable": True,
}

AG_GRID_DEFAULT_OPTIONS = {"pagination": True, "paginationAutoPageSize": True}
