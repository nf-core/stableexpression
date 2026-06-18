process CROSS_JOIN {

    tag "${meta.section} :: ${meta.index_1} vs ${meta.index_2}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

    input:
    tuple val(meta), path("count_chunk_file_1"), path("count_chunk_file_2")

    output:
    tuple val(meta), path('cross_join.*.parquet'),                                                                    emit: data
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    """
    make_cross_join.py \\
        --file1 count_chunk_file_1 \\
        --file2 count_chunk_file_2 \\
        --index1 ${meta.index_1} \\
        --index2 ${meta.index_2}
    """

}
