import plotly.graph_objects as go
import numpy as np
from scipy.stats import gaussian_kde
from dash_extensions.enrich import Input, Output, State, callback

from src.utils.data_management import DataManager

data_manager = DataManager()


##############################################
##############################################
# CALLBACKS
##############################################
##############################################


def register_callbacks():
    @callback(
        Output("sample-counts", "data"),
        Input("sample-dropdown", "value"),
        State("sample-counts", "data"),
        prevent_initial_call=True,
    )
    def update_stored_data(
        sample_dropdown_values: list[str], stored_sample_counts: dict
    ):
        updated_stored_sample_counts = dict(stored_sample_counts)  # deep copy

        # deleting stored data for samples not anymore in the selected list
        for stored_sample in (
            stored_sample_counts
        ):  # we need to copy the list of keys before changing the dict
            if stored_sample not in sample_dropdown_values:
                del updated_stored_sample_counts[stored_sample]

        # storing data for new samples in the selected list
        for sample in sample_dropdown_values:
            if sample not in updated_stored_sample_counts:
                sample_data = data_manager.get_sample_counts(sample)
                updated_stored_sample_counts[sample] = {
                    "counts": sample_data.to_list(),
                    "genes": sample_data.index.to_list(),
                }

        return updated_stored_sample_counts

    @callback(
        Output("sample-graph", "figure"),
        Output("sample-graph", "style"),
        Output("sample_stats_display_accordion_control", "disabled"),
        Output("sample_points_display_accordion_control", "disabled"),
        Input("sample-counts", "data"),
        Input("curve-type", "value"),
        Input("sample-graph-jitter", "value"),
        Input("sample-graph-pointpos", "value"),
        Input("sample-graph-boxmean", "value"),
        Input("sample-graph-display-points", "value"),
        State("sample-graph", "style"),
        prevent_initial_call=True,
    )
    def update_sample_histogram(
        sample_counts: dict,
        curve_type: str,
        jitter: float,
        pointpos: float,
        boxmean: str | bool,
        point_display_mode: str,
        graph_style: dict,
    ):
        if not sample_counts:
            graph_style["display"] = "none"
            return {}, graph_style

        graph_style["display"] = "block"

        fig = go.Figure()

        sample_stats_display_ac_disabled = True
        sample_points_display_ac_disabled = True

        # we need to use the reversed order, otherwise the last traced added is at the top of the graph
        for sample, sample_data in reversed(sample_counts.items()):
            counts = sample_data["counts"]

            if curve_type == "histogram":
                fig.add_trace(go.Histogram(name=sample, x=counts))

            elif curve_type == "kde":
                kde_function = gaussian_kde(counts)
                xvals = np.linspace(min(counts), max(counts), 1000)
                yvals = kde_function(xvals)
                fig.add_trace(go.Scatter(name=sample, x=xvals, y=yvals))

            else:  # boxplot
                # we need to use the reversed order, otherwise the last traced added is at the top of the graph
                fig.add_trace(
                    go.Box(
                        name=sample,
                        x=counts,
                        jitter=jitter,
                        pointpos=pointpos,
                        boxpoints=point_display_mode,
                        boxmean=boxmean,
                        customdata=sample_data["genes"],
                        hovertemplate="Gene: %{customdata}<br>Count: %{x}<br><extra></extra>",
                    )
                )
                # update the layout to remove y-axis labels
                fig.update_layout(yaxis=dict(showticklabels=False))

                sample_stats_display_ac_disabled = False
                sample_points_display_ac_disabled = False

        fig.update_xaxes(range=[0, 1])

        return (
            fig,
            graph_style,
            sample_stats_display_ac_disabled,
            sample_points_display_ac_disabled,
        )
