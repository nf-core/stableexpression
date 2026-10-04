process SEPARATE_COUNTS {

    label 'process_high'

    tag "${meta.platform}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5eb4245644f836b386e2a49a6fbc6f0e6bd810046a54342efce2d2ccd521747/data':
        '' }"

    input:
    tuple val(meta), path(count_file), path(design_file)

    output:
    tuple val(meta), path('*.parquet'), emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    """
    separate_counts_by_batch.py \\
        --counts $count_file \\
        --design $design_file
    """

    stub:
    """
    touch stub.parquet
    """

}
