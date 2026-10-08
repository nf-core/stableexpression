process GLOBAL_STABILITY_SCORE {

    tag "${meta.section}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

    input:
    tuple val(meta), path(platform_stats_score_files)
    path nb_samples_per_platform_file
    val std_penalty_weight
    val null_penalty_weight

    output:
    path "*.stats_with_scores.csv",                                                                                   emit: stats_with_stability_scores
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    """
    compute_global_stability_score.py \\
        --platform-stats-scores "${platform_stats_score_files.join(' ')}" \\
        --nb-samples-per-platform $nb_samples_per_platform_file \\
        --std-penalty-weight $std_penalty_weight \\
        --null-penalty-weight $null_penalty_weight

    mv stats_with_scores.csv ${meta.section}.stats_with_scores.csv
    """

    stub:
    """
    touch section_stub.stats_with_scores.csv
    """

}
