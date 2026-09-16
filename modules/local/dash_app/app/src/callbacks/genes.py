import plotly.graph_objects as go
from dash_extensions.enrich import Input, Output, Serverside, State, callback, ctx
from src.utils.data_management import DataManager

data_manager = DataManager()


##############################################
##############################################
# CALLBACKS
##############################################
##############################################


def get_selected_rows(selected_genes: list[str]) -> list[dict]:
    return data_manager.all_genes_stat_df.filter(
        data_manager.all_genes_stat_df["gene_id"].is_in(selected_genes)
    ).to_dicts()


def register_callbacks():
    @callback(
        Output("gene-counts", "data"),
        Output("gene-dropdown", "value"),
        Output("gene-stats-table", "selectedRows"),
        Input("gene-dropdown", "value"),
        Input("gene-stats-table", "selectedRows"),
        State("gene-counts", "data"),
        # prevent_initial_call=True,
    )
    def update_gene_stored_data(
        selected_genes: list[str], table_selected_rows: list[dict], stored_data: dict
    ) -> dict:
        if ctx.triggered_id == "gene-stats-table":
            # updating selected genes
            if table_selected_rows is not None:
                selected_genes = [row["gene_id"] for row in table_selected_rows]
            else:
                selected_genes = []
        else:
            # ctx.triggered_id is None (callback triggered at app launch / refresh)
            # or ctx.triggered_id == "gene-dropdown":
            # taking the dropdown values as reference (since there is persistence on it)
            table_selected_rows = get_selected_rows(selected_genes)

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

        return Serverside(stored_data), selected_genes, table_selected_rows

    @callback(
        Output("gene-graph", "figure"),
        Output("gene-graph", "style"),
        Input("gene-counts", "data"),
        Input("gene-graph-jitter", "value"),
        Input("gene-graph-pointpos", "value"),
        Input("gene-graph-boxmean", "value"),
        Input("gene-graph-display-points", "value"),
        State("gene-graph", "style"),
        # prevent_initial_call=True,
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
