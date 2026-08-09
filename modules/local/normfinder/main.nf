process NORMFINDER   {

    tag "${meta.section}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/0526f3dbdd23175430f0af81763c1079b3b1425b2cdb2491ab54bb9c0d93d480/data':
        'community.wave.seqera.io/library/numba_numpy_polars_python_tqdm:f42e9bc9f30a29ff' }"

    input:
    tuple val(meta), path(count_file)
    path design_file

    output:
    tuple val(meta), path('stability_values.normfinder.csv'),                                                   emit: stability_values
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions
    tuple val("${task.process}"), val('tqdm'),   eval('python3 -c "import tqdm; print(tqdm.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('numpy'),  eval('python3 -c "import numpy; print(numpy.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('numba'),  eval('python3 -c "import numba; print(numba.__version__)"'),   topic: versions

    script:
    """
    normfinder.py \\
        --counts $count_file \\
        --design $design_file
    """

    stub:

    """
    touch stability_values.normfinder.csv
    """

}
