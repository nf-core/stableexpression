process COMPUTE_STABILITY_SCORES {

    tag "${meta.section}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    tuple val(meta), path(normfinder_stability_file), path(genorm_stability_file)
    path stat_file
    val stability_score_weights

    output:
    path "${meta.section}.stats_with_scores.csv",                                                                     emit: stats_with_stability_scores
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def genorm_stability_file_arg = genorm_stability_file ? "--genorm-stability $genorm_stability_file" : ""
    """
    # limiting number of threads to polars / python
    export POLARS_MAX_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    compute_stability_scores.py \\
        --stats $stat_file \\
        --weights "$stability_score_weights" \\
        --normfinder-stability $normfinder_stability_file \\
        $genorm_stability_file_arg

    mv stats_with_scores.csv ${meta.section}.stats_with_scores.csv
    """

}
