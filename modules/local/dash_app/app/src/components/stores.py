from dash_extensions.enrich import dcc

selected_samples = dcc.Store("selected-sample", storage_type="session")
gene_counts = dcc.Store(id="gene-counts", storage_type="session", data={})
sample_counts = dcc.Store(id="sample-counts", storage_type="session", data={})
# filtered_sample_counts = dcc.Store(id='filtered-sample-counts', storage_type='memory', data={})


stores_to_load = [
    # selected_samples,
    gene_counts,
    sample_counts,
    # filtered_sample_counts
]
