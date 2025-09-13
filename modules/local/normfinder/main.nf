process NORMFINDER   {

    label 'process_single'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0e/0e0445114887dd260f1632afe116b1e81e02e1acc74a86adca55099469b490d9/data':
        'community.wave.seqera.io/library/numba_numpy_polars_tqdm:6923cfab6fc04dec' }"

    input:
    path count_file
    path design_file

    output:
    path('stabilities.csv'),                                                                                            emit: stabilities
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                       topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),       topic: versions

    script:
    """
    normfinder.py \
        --counts $count_file \
        --design $design_file
    """

    stub:
    """
    touch stabilities.csv
    """

}
