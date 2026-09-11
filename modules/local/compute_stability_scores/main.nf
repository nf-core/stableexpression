process COMPUTE_STABILITY_SCORES {

    tag "${meta.section}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

    input:
    tuple val(meta), path(normfinder_stability_file), path(genorm_stability_file), path(section_stat_file)
    val stability_score_weights

    output:
    path "${meta.section}.stats_with_scores.csv",                                                                     emit: stats_with_stability_scores
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def genorm_stability_file_arg = genorm_stability_file ? "--genorm-stability $genorm_stability_file" : ""
    """
    compute_stability_scores.py \\
        --stats $section_stat_file \\
        --weights "$stability_score_weights" \\
        --normfinder-stability $normfinder_stability_file \\
        $genorm_stability_file_arg

    mv stats_with_scores.csv ${meta.section}.stats_with_scores.csv
    """

}
