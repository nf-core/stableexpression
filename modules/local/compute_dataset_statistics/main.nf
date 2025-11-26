process COMPUTE_DATASET_STATISTICS {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5fe497e7a739fa611fedd6f72ab9a3cf925873a5ded3188161fc85fd376b2c1c/data':
        'community.wave.seqera.io/library/pandas_pyarrow_python_scipy:7cad0d297a717147' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta.dataset), path("skewness.txt"),                                                                      topic: skewness
    tuple val(meta.dataset), path("ratio_zeros.txt"),                                                                   topic: ratio_zeros
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                       topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),       topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.dataset}"
    """
    get_dataset_statistics.py \
        --counts $count_file
    """

}
