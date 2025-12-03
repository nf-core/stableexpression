process CLEAN_COUNT_DATA {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    tuple val(meta), path(count_file), path(ks_stats_file)
    val ks_pvalue_threshold

    output:
    tuple val(meta), path('cleaned_counts_filtered.parquet'), optional: true,                                     emit: counts
    tuple val(meta.dataset), path("failure_reason.txt"),      optional: true,                                     topic: clean_count_failure_reason
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    clean_count_data.py \\
        --counts $count_file \\
        --ks-stats $ks_stats_file \\
        --ks-pvalue-threshold $ks_pvalue_threshold
    """

}
