process QUANTILE_NORMALISATION {

    label 'process_low'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/5751f4c7c1eb17d92c2863dec2b7505295e56eafb65ea5a9df66876fbffd24e3/data':
        'community.wave.seqera.io/library/polars_python_scikit-learn:041254a8f0633213' }"

    input:
    tuple val(meta), path(count_file)
    val target_distribution

    output:
    tuple val(meta), path('*.quant_norm.parquet'),                                                                      emit: counts
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                   topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval('python3 -c "import sklearn; print(sklearn.__version__)"'), topic: versions

    script:
    """
    quantile_normalise.py \\
        --counts $count_file \\
        --target-distrib $target_distribution
    """

    stub:
    """
    touch count.cpm.quant_norm.parquet
    """

}
