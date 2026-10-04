process IMPUTE_MISSING_VALUES {

    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/5751f4c7c1eb17d92c2863dec2b7505295e56eafb65ea5a9df66876fbffd24e3/data':
        'community.wave.seqera.io/library/polars_python_scikit-learn:041254a8f0633213' }"

    input:
    path(count_file)

    output:
    tuple val(meta), path('*.corrected.parquet'),                                                                       emit: counts
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                   topic: versions
    tuple val("${task.process}"), val('pyarrow'),      eval('python3 -c "import pyarrow; print(pyarrow.__version__)"'), topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval('python3 -c "import sklearn; print(sklearn.__version__)"'), topic: versions
    tuple val("${task.process}"), val('tqdm'),         eval('python3 -c "import tqdm; print(tqdm.__version__)"'),       topic: versions
    tuple val("${task.process}"), val('joblib'),       eval('python3 -c "import joblib; print(joblib.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('pandas'),       eval('python3 -c "import pandas; print(pandas.__version__)"'),   topic: versions

    script:
    """
    correct_batch_effects_with_recombat.py \\
        --counts $count_file
    """

    stub:
    """
    touch stub.imputed.parquet
    """

}
