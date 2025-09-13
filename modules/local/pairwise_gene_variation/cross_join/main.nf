process CROSS_JOIN {

    // label 'process_medium'
    publishDir "${params.outdir}/pairwise_gene_variation/cross_joins"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    tuple val(meta), path("count_chunk_file_1"), path("count_chunk_file_2")

    output:
    path 'cross_join.*.parquet',                                                                                      emit: data
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    def args = "--task-attempts ${task.attempt}"
    """
    make_cross_join.py --file1 count_chunk_file_1 --file2 count_chunk_file_2 --index1 ${meta.index_1} --index2 ${meta.index_2} $args
    """

}
