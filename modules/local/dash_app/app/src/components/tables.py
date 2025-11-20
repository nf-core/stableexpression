import dash_ag_grid as dag
from src.utils import style
from src.utils.data_management import DataManager

data_manager = DataManager()

NB_GENES_SELECTED_DEFAULT = 10

row_data = data_manager.all_genes_stat_df.to_dicts()
default_selected_rows = data_manager.all_genes_stat_df.head(
    NB_GENES_SELECTED_DEFAULT
).to_dicts()
column_defs = [
    {"field": col, "headerName": col.replace("_", " ").capitalize()}
    for col in data_manager.all_genes_stat_df.columns
]


all_genes_stats_table = dag.AgGrid(
    rowData=row_data,
    columnDefs=column_defs,
    className="ag-theme-alpine",
    # columnSizeOptions=dict(skipHeader=False),
    # columnSize="autoSizetoFit",
    defaultColDef=dict(
        filter=True,
        resizable=True,
        editable=False,
        sortable=True,
    ),
    dashGridOptions=dict(
        pagination=True,
        paginationAutoPageSize=True,
        enableCellTextSelection=True,
        ensureDomOrder=True,
        animateRows=False,
        rowSelection=dict(mode="multiRow"),
        headerCheckboxSelection=False,
        getRowId="params.data.gene_id",
    ),
    selectedRows=default_selected_rows,
    style=style.AG_GRID,
    persistence=True,
    persistence_type="session",
    persisted_props=["selectedRows"],
    id="gene-stats-table",
)
