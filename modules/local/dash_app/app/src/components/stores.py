from dash_extensions.enrich import dcc

selected_samples = dcc.Store("selected-sample", storage_type="session")
gene_counts = dcc.Store(id="gene-counts", storage_type="session", data={})
sample_counts = dcc.Store(id="sample-counts", storage_type="session", data={})
stores_to_load = [
    gene_counts,
    sample_counts,
]
