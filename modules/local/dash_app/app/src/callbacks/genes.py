from dash_extensions.enrich import Input, Output, State, callback
import plotly.graph_objects as go

from src.utils.data_management import DataManager

data_manager = DataManager()


##############################################
##############################################
# CALLBACKS
##############################################
##############################################


def register_callbacks():
    @callback(
        Output("gene-counts", "data"),
        Input("gene-dropdown", "value"),
        State("gene-counts", "data"),
        prevent_initial_call=True,
    )
    def update_gene_stored_data(selected_genes: list[str], stored_data: dict) -> dict:
        # deleting stored data for genes not anymore in the selected list
        for stored_gene in list(
            stored_data.keys()
        ):  # we need to copy the list of keys before changing the dict
            if stored_gene not in selected_genes:
                del stored_data[stored_gene]

        # storing data for new genes in the selected list
        for gene in selected_genes:
            if gene not in stored_data:
                gene_data = data_manager.get_gene_counts(gene)
                stored_data[gene] = {
                    "counts": gene_data.to_list(),
                    "samples": gene_data.index.to_list(),
                }
        return stored_data

    @callback(
        Output("gene-graph", "figure"),
        Output("gene-graph", "style"),
        Input("gene-counts", "data"),
        Input("gene-graph-jitter", "value"),
        Input("gene-graph-pointpos", "value"),
        Input("gene-graph-boxmean", "value"),
        Input("gene-graph-display-points", "value"),
        State("gene-graph", "style"),
        prevent_initial_call=True,
    )
    def update_gene_graph(
        gene_stored_data: dict,
        jitter: float,
        pointpos: float,
        boxmean: str | bool,
        point_display_mode: str,
        graph_style: dict,
    ):
        if not gene_stored_data:
            graph_style["display"] = "none"
            return {}, graph_style

        graph_style["display"] = "block"

        fig = go.Figure()

        # we need to use the reversed order, otherwise the last traced added is at the top of the graph
        for gene, gene_data in reversed(gene_stored_data.items()):
            fig.add_trace(
                go.Box(
                    name=gene,
                    x=gene_data["counts"],
                    boxmean=boxmean,
                    jitter=jitter,
                    pointpos=pointpos,
                    boxpoints=point_display_mode,
                    customdata=gene_data["samples"],
                    hovertemplate="Sample: %{customdata}<br>Normalised count: %{x}<br><extra></extra>",
                    showlegend=False,
                )
            )

        fig.update_layout(xaxis=dict(range=[0, 1]), yaxis=dict(ticklabelstandoff=10))

        return fig, graph_style
