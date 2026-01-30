process FILTER_OUT_SAMPLES_WITH_TOO_MANY_ZEROS {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    tuple val(meta), path(count_file)
    val(max_zero_ratio)

    output:
    tuple val(meta), path("*.zeros_filtered.parquet"), optional: true,                                            emit: counts
    tuple val(meta.dataset), env("NB_KEPT_SAMPLES"), env("NB_REJECTED_SAMPLES"),                                  topic: stats_zero_values_filter
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    """
    filter_out_samples_with_too_many_zeros.py \\
        --counts $count_file \\
        --max-zero-ratio $max_zero_ratio

    NB_REJECTED_SAMPLES=\$(cat nb_rejected_samples.csv)
    NB_KEPT_SAMPLES=\$(cat nb_kept_samples.csv)
    """

}
