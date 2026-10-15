process GENORM {

    tag "${meta.platform} :: ${meta.section}"
    label 'process_medium_requirements'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83ffc025ce0d913f9eaa4e786c9ebd2817ddc57f85d6596663f0e8a290872fab/data':
        'community.wave.seqera.io/library/polars_python:e0e89fee0d134a04' }"

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
