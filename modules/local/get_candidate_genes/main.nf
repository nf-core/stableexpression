process GET_CANDIDATE_GENES {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_file
    path stat_file
    val candidate_selection_descriptor
    val nb_most_stable_genes
    val min_pct_quantile_expr_level

    output:
    path 'candidate_counts.parquet',                                                                                  emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    get_candidate_genes.py \\
        --counts $count_file \\
        --stats $stat_file \\
        --candidate_selection_descriptor $candidate_selection_descriptor \\
        --nb-top-stable-genes $nb_most_stable_genes \\
        --min-pct-quantile-expr-level $min_pct_quantile_expr_level
    """

}
