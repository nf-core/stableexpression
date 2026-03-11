import dash_mantine_components as dmc
from dash_iconify import DashIconify
from src.components import graphs, tables
from src.utils import style

gene_icon = DashIconify(icon="material-symbols:genetics", width=20)

sample_icon = DashIconify(icon="ic:baseline-dashboard-customize", width=20)


tabs = dmc.Tabs(
    children=[
        dmc.TabsList(
            children=[
                dmc.TabsTab(
                    dmc.Text("Counts / gene", fw=500),
                    className="genes-tabitem",
                    color="teal",
                    leftSection=gene_icon,
                    value="genes",
                    style=style.HEADER_TABLIST_ITEM,
                ),
                dmc.TabsTab(
                    dmc.Text("Counts / sample", fw=500),
                    className="samples-tabitem",
                    leftSection=sample_icon,
                    value="samples",
                    color="red",
                    style=style.HEADER_TABLIST_ITEM,
                ),
                dmc.TabsTab(
                    dmc.Text("Statistics - all genes", fw=500),
                    leftSection=sample_icon,
                    value="gene_stats",
                    color="orange",
                    style=style.HEADER_TABLIST_ITEM,
                ),
            ],
            style=style.HEADER_TABLIST,
        ),
        dmc.TabsPanel(
            children=[
                dmc.Text("dhkhg"),
                graphs.gene_graph,
            ],
            style=style.TABS_PANEL,
            value="genes",
        ),
        dmc.TabsPanel(
            children=[
                graphs.sample_graph,
            ],
            style=style.TABS_PANEL,
            value="samples",
        ),
        dmc.TabsPanel(
            children=[tables.all_genes_stats_table],
            style=style.TABS_PANEL,
            value="gene_stats",
        ),
    ],
    id="tabs",
    variant="default",
    radius="md",
    orientation="horizontal",
    placement="right",
    value="genes",
    persistence=True,
    persisted_props=["value"],
    persistence_type="session",
    style=style.TAB,
)

settings_button = dmc.Button(
    "Select data / options",
    id="settings-button",
    className="settings-button",
    color="teal",
    style=style.SETTINGS_BUTTON,
)

header = dmc.Grid(
    children=[
        dmc.GridCol(tabs, span=10),
        dmc.GridCol(
            settings_button, span=2, style={"textAlign": "right", "marginTop": "20px"}
        ),
    ],
    style={"marginRight": "20px"},
    # gutter="xl",
)
