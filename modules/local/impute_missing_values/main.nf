process IMPUTE_MISSING_VALUES {

    label 'process_high'
    label 'process_long'

    tag "${meta.platform}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a510353cae7c5f5915ec9c3b3de4af40c8b90e45b1b41ae851581cebba2f563/data':
        'community.wave.seqera.io/library/python_polars_scikit-learn_joblib:5bc5a48bc50fde19' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.imputed.parquet'),                                                                         emit: counts
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                   topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval('python3 -c "import sklearn; print(sklearn.__version__)"'), topic: versions
    tuple val("${task.process}"), val('joblib'),       eval('python3 -c "import joblib; print(joblib.__version__)"'),   topic: versions

    script:
    """
    impute_missing_values.py \\
        --counts $count_file \\
        --out ${meta.platform}.imputed.parquet \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch stub.imputed.parquet
    """

}
