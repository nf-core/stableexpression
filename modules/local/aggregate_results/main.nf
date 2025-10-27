process AGGREGATE_RESULTS {

    label 'process_single'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_file
    path stat_file
    path rnaseq_dataset_stat_file, stageAs: "*/*"
    path microarray_dataset_stat_file, stageAs: "*/*"
    path metadata_files, stageAs: "*/*"
    path mapping_files, stageAs: "*/*"

    output:
    path 'top_stable_genes_summary.csv',                                                                              emit: top_stable_genes_summary
    path 'stats_all_genes.csv',                                                                                       emit: stats_all_genes
    path 'all_counts_filtered.parquet',                                                                               emit: all_counts_filtered
    path 'top_stable_genes_transposed_counts_filtered.csv',                                                           emit: top_stable_genes_transposed_counts_filtered
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def rnaseq_dataset_stat_file_arg = rnaseq_dataset_stat_file ? "--rnaseq $rnaseq_dataset_stat_file" : ""
    def microarray_dataset_stat_file_arg = microarray_dataset_stat_file ? "--microarray $microarray_dataset_stat_file" : ""
    """
    aggregate_results.py \\
        --counts $count_file \\
        --stats $stat_file \\
        --metadata "$metadata_files" \\
        --mappings "$mapping_files" \\
        $rnaseq_dataset_stat_file_arg \\
        $microarray_dataset_stat_file_arg \\
    """

}
