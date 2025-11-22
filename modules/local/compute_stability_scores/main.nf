process COMPUTE_STABILITY_SCORES {

    label 'process_high'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/01/0118e0577564644b18f94fa6525fe3a2aec845721081b55d82a18e803a50ab17/data':
        'community.wave.seqera.io/library/polars_scikit-learn:036e189d7c1f9704' }"

    input:
    path stat_file
    val stability_score_weights
    path normfinder_stability_file
    val genorm_stability_file

    output:
    path 'stats_with_scores.csv',                                                                                     emit: stats_with_stability_scores
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def genorm_stability_file_arg = genorm_stability_file ? "--genorm-stability $genorm_stability_file" : ""
    """
    compute_stability_scores.py \\
        --stats $stat_file \\
        --weights "$stability_score_weights" \\
        --normfinder-stability $normfinder_stability_file \\
        $genorm_stability_file_arg
    """

}
