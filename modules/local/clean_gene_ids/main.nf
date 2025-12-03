process CLEAN_GENE_IDS {

    label 'process_low'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9b43e446f2c3b794644fd4c1c86ab09ba0afafc0c02e3fcdf45509ffc89fc4d/data':
        'community.wave.seqera.io/library/pandas_polars:29ea1468b5490a67' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.cleaned.csv'),                 optional: true,                                           emit: counts
    tuple val(meta.dataset), path("failure_reason.txt"),    optional: true,                                           topic: id_cleaning_failure_reason
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    clean_gene_ids.py \\
        --count-file "$count_file"
    """


    stub:
    """
    touch fake.cleaned.csv
    """

}
