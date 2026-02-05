process IMPUTE_MISSING_VALUES {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb8feda3812519f6f6f085e1d058f534b0aedba570c1443c4479d79975e81906/data':
        'community.wave.seqera.io/library/polars_scikit-learn:a30d22b117dad962' }"

    input:
    tuple val(meta), path(count_file)
    val missing_value_imputer

    output:
    tuple val(meta), path('*.imputed.parquet'),                                                                         emit: counts
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                   topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval('python3 -c "import sklearn; print(sklearn.__version__)"'), topic: versions

    script:
    """
    # limiting number of threads used by polars
    export POLARS_MAX_THREADS=${task.cpus}

    impute_missing_values.py \\
        --counts $count_file \\
        --imputer $missing_value_imputer \\
        --memory "${task.memory}"
    """

}
