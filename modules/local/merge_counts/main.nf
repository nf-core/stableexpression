process MERGE_COUNTS {

    tag "${meta.platform}"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/90617e987f709570820b8e7752baf9004ba85917111425d4b44b429b27b201ca/data':
        'community.wave.seqera.io/library/polars_tqdm:54b124dde91d1bf3' }"

    input:
    tuple val(meta), path(count_files, stageAs: "?/*")

    output:
    tuple val(meta), path('all_counts.parquet'),                                                                      emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('tqdm'),     eval('python3 -c "import tqdm; print(tqdm.__version__)"'),         topic: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    merge_counts.py \\
        --counts "$count_files"
    """

}
