import dash_ag_grid as dag

from src.utils import style


from src.utils.data_management import DataManager

data_manager = DataManager()


def format_col_name(col: str):
    return col.replace("_", " ").capitalize()


def format_col_width(col):
    return col.replace("_", " ").capitalize()


gene_stats_table = dag.AgGrid(
    rowData=data_manager.stat_df.to_dicts(),
    columnDefs=[
        {"field": col, "headerName": format_col_name(col)}
        for col in data_manager.stat_df.columns
    ],
    className="ag-theme-alpine",
    columnSizeOptions=dict(
        skipHeader=False,
    ),
    # columnSizeOptions={'defaultMinWidth': 300},
    columnSize="autoSize",
    defaultColDef=dict(
        # type='rightAligned',
        filter=True,
        resizable=True,
        editable=False,
        sortable=True,
    ),
    dashGridOptions=dict(
        # domLayout="autoHeight",
        pagination=True,
        paginationAutoPageSize=True,
        enableCellTextSelection=True,
        enableRangeSelection=True,
    ),
    style=style.AG_GRID,
    id="gene-stats-table",
)
