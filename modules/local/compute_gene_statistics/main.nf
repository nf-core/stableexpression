process COMPUTE_GENE_STATISTICS {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a1de3eb1a051ef4527661296f6a3165c7d1b0fd8707d21844bfad6483dce4dcb/data':
        'community.wave.seqera.io/library/polars_python:100fa0b0355e4749' }"

    input:
    tuple val(meta), path(count_file), path(imputed_count_file)
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
        args += " --imputed-counts $imputed_count_file"
    }
    """
    compute_gene_statistics.py \\
        --counts $count_file \\
        --ratio-nulls-per-sample $ratio_nulls_per_samples \\
        --max-ratio-null-valid-sample $max_null_ratio_valid_sample \\
        $args
    """

}
