process GENORM {

    tag "${meta.platform} :: ${meta.section}"
    label 'process_medium_requirements'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/0526f3dbdd23175430f0af81763c1079b3b1425b2cdb2491ab54bb9c0d93d480/data':
        'community.wave.seqera.io/library/numba_numpy_polars_python_tqdm:f42e9bc9f30a29ff' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path("*.m_measures.csv"),                                                                  emit: m_measures
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    """
    genorm.py \\
        --file "$count_file" \\
        --out ${meta.section}.m_measures.csv
    """

    stub:
    """
    touch stub.m_measures.csv
    """

}
