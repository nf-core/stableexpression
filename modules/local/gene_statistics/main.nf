process GENE_STATISTICS {

    label 'process_low'
    debug true

    publishDir "${params.outdir}/gene_statistics"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/polars_python:1636dad2fac97d0b':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_file
    path metadata_files, stageAs: "?/*"
    path mapping_files, stageAs: "?/*"
    val m_measure_file

    output:
    path 'top_stable_genes_summary.csv',                                                                              emit: top_stable_genes_summary
    path 'stats_all_genes.csv',                                                                                       emit: all_statistics
    path 'all_log_counts.csv',                                                                                        emit: log_counts
    path 'top_stable_genes_transposed_log_counts.csv',                                                                emit: top_stable_genes_transposed_log_counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    """
    get_gene_statistics.py \
        --counts $count_file \
        --metadata "$metadata_files" \
        --mappings "$mapping_files" \
        --m-measures $m_measure_file
    """

}
