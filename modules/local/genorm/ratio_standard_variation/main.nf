process RATIO_STANDARD_VARIATION {

    tag "${meta.section} :: ${meta.index_1} vs ${meta.index_2}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path('std.*.parquet'),                                                                           emit: data
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    def args = "--task-attempts ${task.attempt}"
    """
    # limiting number of threads to polars / python
    export POLARS_MAX_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    get_ratio_standard_variation.py --file $file $args
    """

}
