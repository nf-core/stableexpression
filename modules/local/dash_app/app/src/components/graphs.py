from dash_extensions.enrich import dcc

from src.utils import style


def get_graph(graph_id: str):
    return dcc.Graph(id=graph_id, figure={}, style=style.GRAPH)


gene_graph = get_graph("gene-graph")

sample_graph = get_graph("sample-graph")
