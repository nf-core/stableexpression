import dash_mantine_components as dmc
from src.utils import style
from src.utils.data_management import DataManager

data_manager = DataManager()

sorted_samples = data_manager.get_sorted_samples()

NB_SAMPLES_DEFAULT = 10

sample_selection_stack = dmc.Stack(
    [
        dmc.MultiSelect(
            id="sample-dropdown",
            label="Select list of samples to visualise",
            placeholder="Select samples",
            nothingFoundMessage="No samples found",
            data=sorted_samples,
            value=sorted_samples[:NB_SAMPLES_DEFAULT],
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
            # style=style.DROPDOWN,
            comboboxProps={
                "shadow": "md",
                "transitionProps": {"transition": "pop", "duration": 200},
            },
        )
    ],
    align="left",
    gap="xl",
)


sample_graph_plot_type_stack = dmc.Stack(
    [
        dmc.Text("Type of plot", style=style.STACK_SUBSECTION_TITLE),
        dmc.SegmentedControl(
            id="curve-type",
            value="ng",
            color="teal",
            data=[
                {"value": "histogram", "label": "Histogram"},
                {"value": "kde", "label": "Kde"},
                {"value": "boxplot", "label": "Box-plot"},
            ],
            mb=10,
        ),
    ],
    align="left",
    gap="xl",
)


sample_graph_stats_display_stack = dmc.Stack(
    [
        dmc.Text(
            "Display mean / standard deviation", style=style.STACK_SUBSECTION_TITLE
        ),
        dmc.SegmentedControl(
            id="sample-graph-boxmean",
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

sample_graph_points_display_stack = dmc.Stack(
    [
        dmc.Text("Display points", style=style.STACK_SUBSECTION_TITLE),
        dmc.SegmentedControl(
            id="sample-graph-display-points",
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
            id="sample-graph-pointpos",
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
            id="sample-graph-jitter",
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
    value="sample_selection",
    children=[
        dmc.AccordionItem(
            [
                dmc.AccordionControl("Sample selection"),
                dmc.AccordionPanel(sample_selection_stack),
            ],
            value="sample_selection",
        ),
        dmc.AccordionItem(
            [
                dmc.AccordionControl(
                    "Plot customisation",
                    id="sample_plot_customisation_accordion_control",
                ),
                dmc.AccordionPanel(sample_graph_plot_type_stack),
            ],
            value="sample_plot_customisation",
        ),
        dmc.AccordionItem(
            [
                dmc.AccordionControl(
                    "Statistics display", id="sample_stats_display_accordion_control"
                ),
                dmc.AccordionPanel(sample_graph_stats_display_stack),
            ],
            value="sample_stats_display",
        ),
        dmc.AccordionItem(
            [
                dmc.AccordionControl(
                    "Points display", id="sample_points_display_accordion_control"
                ),
                dmc.AccordionPanel(sample_graph_points_display_stack),
            ],
            value="sample_points_display",
        ),
    ],
    id="sidebar-samples-items",
    style={"marginTop": "20px", "display": "none"},
)
