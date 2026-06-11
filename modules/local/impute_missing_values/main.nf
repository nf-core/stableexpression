process IMPUTE_MISSING_VALUES {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/5751f4c7c1eb17d92c2863dec2b7505295e56eafb65ea5a9df66876fbffd24e3/data':
        'community.wave.seqera.io/library/polars_python_scikit-learn:041254a8f0633213' }"

    input:
    tuple val(meta), path(count_file)
    val missing_value_imputer
    val knn_n_neighbours
    val iterative_max_iter
    val iterative_n_nearest_features

    output:
    tuple val(meta), path('*.imputed.parquet'),                                                                         emit: counts
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                   topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval('python3 -c "import sklearn; print(sklearn.__version__)"'), topic: versions

    script:
    """
    impute_missing_values.py \\
        --counts $count_file \\
        --imputer $missing_value_imputer \\
        --knn-n-neighbours $knn_n_neighbours \\
        --iterative-max-iter $iterative_max_iter \\
        --iterative-n-nearest-features $iterative_n_nearest_features
    """

}
