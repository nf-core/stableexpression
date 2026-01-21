process NORMFINDER   {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0e/0e0445114887dd260f1632afe116b1e81e02e1acc74a86adca55099469b490d9/data':
        'community.wave.seqera.io/library/numba_numpy_polars_tqdm:6923cfab6fc04dec' }"

    input:
    path count_file
    path design_file

    output:
    path('stability_values.normfinder.csv'),                                                                            emit: stability_values
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                       topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),       topic: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    normfinder.py \
        --counts $count_file \
        --design $design_file
    """

    stub:

    """
    touch stability_values.normfinder.csv
    """

}
