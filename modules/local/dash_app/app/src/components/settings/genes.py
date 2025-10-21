import dash_mantine_components as dmc

from src.utils import style
from src.utils.data_management import DataManager

data_manager = DataManager()

gene_selection_stack = dmc.Stack(
    [
        dmc.MultiSelect(
            id="gene-dropdown",
            label=dmc.Text("Genes to display", fw=600, style={"paddingBottom": "5px"}),
            placeholder="Select genes of interest",
            nothingFoundMessage="No gene found",
            data=data_manager.genes,
            value=[],
            w=400,
            clearable=True,
            searchable=True,
            limit=100,
            maxValues=20,
            size="sm",
            checkIconPosition="right",
            hidePickedOptions=True,
            disabled=False,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            style=style.DROPDOWN,
            comboboxProps={
                "shadow": "md",
                "transitionProps": {"transition": "pop", "duration": 200},
            },
        )
    ],
    align="left",
    gap="xl",
)

gene_graph_stats_display_stack = dmc.Stack(
    [
        dmc.Text(
            "Display mean / standard deviation", style=style.STACK_SUBSECTION_TITLE
        ),
        dmc.SegmentedControl(
            id="gene-graph-boxmean",
            value="sd",
            color="teal",
            data=[
                {"value": False, "label": "None"},
                {"value": True, "label": "Mean only"},
                {"value": "sd", "label": "Mean + Std"},
            ],
            mb=10,
        ),
    ],
    align="left",
    gap="xl",
)

gene_graph_points_display_stack = dmc.Stack(
    [
        dmc.Text("Display points", style=style.STACK_SUBSECTION_TITLE),
        dmc.SegmentedControl(
            id="gene-graph-display-points",
            value="outliers",
            color="teal",
            data=[
                {"value": "outliers", "label": "Outliers"},
                {"value": "suspectedoutliers", "label": "Suspected Outliers"},
                {"value": "all", "label": "All points"},
            ],
            mb=10,
        ),
        dmc.Text(
            "Position of points relatively to boxes", style=style.STACK_SUBSECTION_TITLE
        ),
        dmc.Slider(
            id="gene-graph-pointpos",
            value=-1.8,
            color="teal",
            min=-2,
            max=2,
            step=0.1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
        dmc.Text(
            "Spreading of displayed points (jitter)", style=style.STACK_SUBSECTION_TITLE
        ),
        dmc.Slider(
            id="gene-graph-jitter",
            value=0.3,
            color="teal",
            min=0,
            max=1,
            step=0.1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
    ],
    align="left",
    gap="xl",
)

sidebar_stack = dmc.Accordion(
    value="gene_selection",
    children=[
        dmc.AccordionItem(
            [
                dmc.AccordionControl("Gene selection"),
                dmc.AccordionPanel(gene_selection_stack),
            ],
            value="gene_selection",
        ),
        dmc.AccordionItem(
            [
                dmc.AccordionControl("Statistics display"),
                dmc.AccordionPanel(gene_graph_stats_display_stack),
            ],
            value="gene_stats_display",
        ),
        dmc.AccordionItem(
            [
                dmc.AccordionControl("Points display"),
                dmc.AccordionPanel(gene_graph_points_display_stack),
            ],
            value="gene_points_display",
        ),
    ],
    id="sidebar-genes-items",
    style={"marginTop": "20px", "display": "none"},
)
