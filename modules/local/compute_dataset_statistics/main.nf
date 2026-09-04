process COMPUTE_DATASET_STATISTICS {

    label 'process_single_cpu'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/878943dcc1b8e30cd535a41886e0f75fcd8bbe667b2d2b0bc4adb0c549539e64/data':
        'community.wave.seqera.io/library/polars_python:07cce0ec1b0aeb84' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta.dataset), path("skewness.txt"),                                                                      topic: skewness
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                       topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),       topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.dataset}"
    """
    compute_dataset_statistics.py \\
        --counts $count_file
    """

}
