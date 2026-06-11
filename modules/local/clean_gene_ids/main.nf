process CLEAN_GENE_IDS {

    label 'process_low'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/878943dcc1b8e30cd535a41886e0f75fcd8bbe667b2d2b0bc4adb0c549539e64/data':
        'community.wave.seqera.io/library/polars_python:07cce0ec1b0aeb84' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.cleaned.parquet'),             optional: true,                                           emit: counts
    tuple val(meta.dataset), path("failure_reason.txt"),    optional: true,                                           topic: id_cleaning_failure_reason
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    """
    clean_gene_ids.py \\
        --count-file "$count_file"
    """

}
