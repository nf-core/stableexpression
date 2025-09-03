process COMPUTE_GLOBAL_GENE_STATISTICS {

    label 'process_low'

    errorStrategy = {
        if (task.exitStatus == 100) {
            log.error(
                "No count could be found before merging datasets! "
                + "Please check the provided accessions and datasets and run again"
                )
            return 'terminate'
        }
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_file
    path platform_statistic_files, stageAs: "?/*"
    path metadata_files, stageAs: "?/*"
    path mapping_files, stageAs: "?/*"
    val nb_top_stable_genes

    output:
    path 'top_stable_genes_summary.csv',                                                                              emit: top_stable_genes_summary
    path 'stats_all_genes.csv',                                                                                       emit: all_statistics
    path 'top_stable_genes_transposed_counts_filtered.csv',                                                           emit: top_stable_genes_transposed_counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    compute_final_gene_statistics.py \\
        --counts $count_file \\
        --stats "$platform_statistic_files" \\
        --metadata "$metadata_files" \\
        --mappings "$mapping_files" \\
        --nb-top-stable-genes $nb_top_stable_genes
    """

}
