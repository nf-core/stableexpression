process QUANTILE_NORMALISATION {

    label 'process_low'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2df931a4ea181fe1ea9527abe0fd4aff9453d6ea56d56aee7c4ac5dceed611e3/data':
        'community.wave.seqera.io/library/pandas_pyarrow_python_scikit-learn:6f85e3c4d1706e81' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.quant_norm.parquet'),                                                                      emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                       topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),       topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval('python3 -c "import sklearn; print(sklearn.__version__)"'), topic: versions
    tuple val("${task.process}"), val('pyarrow'), eval('python3 -c "import pyarrow; print(pyarrow.__version__)"'),      topic: versions


    script:
    """
    quantile_normalise.py --counts $count_file
    """


    stub:
    """
    touch count.cpm.quant_norm.parquet
    """

}
