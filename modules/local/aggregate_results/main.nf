process AGGREGATE_RESULTS {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_file
    path stat_file
    path platform_stat_files, stageAs: "?/*"
    path metadata_files
    path mapping_files

    output:
    path 'all_genes_summary.csv',                                                                                     emit: all_genes_summary
    path 'most_stable_genes_summary.csv',                                                                             emit: most_stable_genes_summary
    path 'all_counts_filtered.parquet',                                                                               emit: all_counts_filtered
    path 'most_stable_genes_transposed_counts_filtered.csv',                                                          emit: most_stable_genes_transposed_counts_filtered
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def mapping_files_arg = mapping_files ? "--mappings " + "$mapping_files" : ""
    def metadata_files_arg = metadata_files ? "--metadata " + "$metadata_files" : ""
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    aggregate_results.py \\
        --counts $count_file \\
        --stats $stat_file \\
        --platform-stats $platform_stat_files \\
        $mapping_files_arg \\
        $metadata_files_arg
    """

}
