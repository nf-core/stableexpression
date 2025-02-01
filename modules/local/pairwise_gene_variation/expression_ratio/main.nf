process EXPRESSION_RATIO {

    // label 'process_medium'
    publishDir "${params.outdir}/pairwise_gene_variation/expression_ratios"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars_python:1636dad2fac97d0b':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path file

    output:
    path 'ratios.*.parquet',                                                                                            emit: data
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    def args = "--task-attempts ${task.attempt}"
    """
    make_pairwise_gene_expression_ratio.py --file $file
    """

}
