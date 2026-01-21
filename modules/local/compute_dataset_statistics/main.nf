process COMPUTE_DATASET_STATISTICS {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta.dataset), path("skewness.txt"),                                                                      topic: skewness
    tuple val(meta.dataset), path("ratio_zeros.txt"),                                                                   topic: ratio_zeros
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                       topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),       topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.dataset}"
    """
    compute_dataset_statistics.py \\
        --counts $count_file
    """

}
