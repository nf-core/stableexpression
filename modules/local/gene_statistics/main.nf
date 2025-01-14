process GENE_STATISTICS {

    label 'process_low'

    publishDir "${params.outdir}/gene_statistics"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pandas_polars_python:2130fd70334a9168':
        'community.wave.seqera.io/library/pandas_polars_python:d34179f2767e1bc7' }"

    input:
    path count_files, stageAs: "?/*"
    path metadata_files, stageAs: "?/*"
    path mapping_files, stageAs: "?/*"

    output:
    path 'stats_all_genes.csv',                                                                                       emit: stats_all_genes
    path 'stats_most_stable_genes.csv',                                                                               emit: stats_most_stable_genes
    path 'log_count_summary.csv',                                                                                     emit: log_count_summary
    path 'count_summary.parquet',                                                                                     emit: count_summary
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    """
    merge_counts_and_get_gene_statistics.py --counts "$count_files" --metadata "$metadata_files" --mappings "$mapping_files"
    """

}
