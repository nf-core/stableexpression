process COMPUTE_GENE_STATISTICS {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

    input:
    tuple val(meta), path(count_file, name: 'count_file.parquet'), path(imputed_count_file, name: 'imputed_count_file.parquet')
    path ratio_nulls_per_samples
    val max_null_ratio_valid_sample

    output:
    path '*stats_all_genes.csv',                                                                                  emit: stats
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    def args = task.ext.args ?: ''
    if ( meta.platform != "all" ) {
        args += " --platform $meta.platform"
    }
    if ( imputed_count_file ) {
        args += " --imputed-counts imputed_count_file.parquet"
    }
    """
    compute_gene_statistics.py \\
        --counts count_file.parquet \\
        --ratio-nulls-per-sample $ratio_nulls_per_samples \\
        --max-ratio-null-valid-sample $max_null_ratio_valid_sample \\
        $args
    """

}
