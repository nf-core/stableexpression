process FILTER_GENES {

    label 'process_high'

    tag "${meta.platform}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

    input:
    tuple val(meta), path(count_file)
    val max_zero_ratio

    output:
    tuple val(meta), path("*.genes_filtered.parquet"),                                                          emit: counts
    tuple val(meta.platform), env("NB_KEPT_GENES"), env("NB_REJECTED_GENES"),                                   topic: mqc_filtered_genes
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    """
    filter_genes.py \\
        --counts $count_file \\
        --max-zero-ratio $max_zero_ratio

    NB_KEPT_GENES=\$(cat nb_kept_genes.csv)
    NB_REJECTED_GENES=\$(cat nb_rejected_genes.csv)
    """

    stub:
    """
    touch stub.genes_filtered.parquet
    NB_KEPT_GENES=0
    NB_REJECTED_GENES=0
    """

}
