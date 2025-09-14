process COMPUTE_STABILITY_SCORES {

    label 'process_single'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path stat_file
    path stability_files, stageAs: "?/*"
    val scoring_base

    output:
    path 'stats_with_scores.csv',                                                                                     emit: stats_with_stability_scores
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    """
    compute_stability_scores.py \\
        --stats $stat_file \\
        --stabilities "$stability_files" \\
        --scoring-base $scoring_base
    """

}
