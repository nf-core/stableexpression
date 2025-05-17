process DATASET_STATISTICS {

    label 'process_low'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5fe497e7a739fa611fedd6f72ab9a3cf925873a5ded3188161fc85fd376b2c1c/data':
        'community.wave.seqera.io/library/pandas_pyarrow_python_scipy:7cad0d297a717147' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.dataset_stats.csv'),                                                                       emit: stats
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                       topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),       topic: versions
    tuple val("${task.process}"), val('scipy'),    eval('python3 -c "import scipy; print(scipy.__version__)"'),         topic: versions
    tuple val("${task.process}"), val('pyarrow'),  eval('python3 -c "import pyarrow; print(pyarrow.__version__)"'),     topic: versions


    script:
    """
    get_dataset_statistics.py --counts $count_file
    """


    stub:
    """
    touch count.cpm.dataset_stats.csv
    """

}
