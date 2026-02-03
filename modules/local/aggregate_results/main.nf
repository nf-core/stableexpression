process AGGREGATE_RESULTS {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a040ba30cbb433a3a6e84ca9881dce77e23762a2b860bdea21b252296a366d20/data':
        'community.wave.seqera.io/library/polars_python_pyyaml:8b53cd142171d9f8' }"

    input:
    path count_file
    path stat_score_files
    path platform_stat_files, stageAs: "?/*"
    path metadata_files
    path mapping_files
    path multiqc_config

    output:
    path 'all_genes_summary.csv',                                                                                     emit: all_genes_summary
    path '*most_stable_genes_summary.csv',                                                                            emit: most_stable_genes_summary
    path '*most_stable_genes_transposed_counts.csv',                                                                  emit: most_stable_genes_transposed_counts_filtered
    path 'custom_content_multiqc_config.yaml',                                                                        emit: custom_content_multiqc_config
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('pyyaml'),   eval('python3 -c "import yaml; print(yaml.__version__)"'),         topic: versions

    script:
    def mapping_files_arg = mapping_files   ? "--mappings " + "$mapping_files"  : ""
    def metadata_files_arg = metadata_files ? "--metadata " + "$metadata_files" : ""
    """
    aggregate_results.py \\
        --cpus ${task.cpus} \\
        --memory "${task.memory}" \\
        --counts $count_file \\
        --stats-with-scores $stat_score_files \\
        --platform-stats $platform_stat_files \\
        --multiqc-config $multiqc_config \\
        $mapping_files_arg \\
        $metadata_files_arg
    """

}
