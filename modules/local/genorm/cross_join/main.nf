process CROSS_JOIN {

    tag "${meta.platform} :: ${meta.section}"
    label 'process_low_requirements'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path('*.cross_join.parquet'), emit: data
    tuple val("GENORM_CROSS_JOIN"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("GENORM_CROSS_JOIN"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    """
    make_cross_join.py \\
        --file $counts
    """

    stub:
    """
    touch stub.cross_join.parquet
    """

}
