process COMPUTE_DATASET_STATISTICS {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

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
