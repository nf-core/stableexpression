process FILTER_OUT_SAMPLES_WITH_TOO_MANY_MISSING_VALUES {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

    input:
    tuple val(meta), path(count_file)
    path valid_gene_ids
    val max_null_ratio

    output:
    tuple val(meta), path("*.nulls_filtered.parquet"), optional: true, emit: counts
    path("ratio_null_values_per_sample.csv"),                          emit: ratio_nulls_per_sample
    tuple val(meta.dataset), path("ratio_null_values.csv"),                                                                                     topic: ratio_nulls
    tuple val(meta.dataset), env("NB_KEPT_SAMPLES"), env("NB_REJECTED_SAMPLES"),                                                                topic: mqc_missing_values_filter_stats
    tuple val("FILTER_OUT_SAMPLES_WITH_TOO_MANY_MISSING_VALUES"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("FILTER_OUT_SAMPLES_WITH_TOO_MANY_MISSING_VALUES"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    """
    filter_out_samples_with_too_many_missing_values.py \\
        --counts $count_file \\
        --valid-gene-ids $valid_gene_ids \\
        --max-null-ratio $max_null_ratio

    NB_REJECTED_SAMPLES=\$(cat nb_rejected_samples.csv)
    NB_KEPT_SAMPLES=\$(cat nb_kept_samples.csv)
    """

}
