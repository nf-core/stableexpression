process COMPUTE_GENE_STATISTICS_PER_PLATFORM {

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
    val platform

    output:
    path '*_stats_all_genes.csv',                                                                                       emit: stats
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    """
    compute_gene_statistics_per_platform.py \\
        --counts $count_file \\
        --platform $platform
    """

}
