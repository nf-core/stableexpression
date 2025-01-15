process PAIRWISE_GENE_VARIATION {

    // label 'process_medium'

    publishDir "${params.outdir}/pairwise_gene_variation"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars_psutil_python_tqdm:ebc5651c8258f46d':
        'community.wave.seqera.io/library/polars_psutil_python_tqdm:e7cefc1d3226eea4' }"

    input:
    path count_file

    output:
    path 'm_measures.csv',                                                                                            emit: m_measures
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('psutil'),   eval('python3 -c "import psutil; print(psutil.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('tqdm'),     eval('python3 -c "import tqdm; print(tqdm.__version__)"'),         topic: versions


    script:
    def args = "--task-attempts ${task.attempt}"
    """
    run_pairwise_gene_variation_analysis.py --counts "$count_file" $args
    """

}
