process CLEAN_COUNT_DATA {

    label 'process_single'

    errorStrategy {
        if (task.exitStatus == 101) {
            /*
            log.warning(
                "No more valid sample after checking p-value of Kolmogorow-Smirnoff test against target distribution! "
                + "You can try a more flexible approach by setting again the value of the ks_pvalue_threshold parameter. "
                + "Provide a negative value to disable this filter."
            )
            */
            return 'ignore'
        }
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    tuple val(meta), path(count_file), path(ks_stats_file)
    val ks_pvalue_threshold

    output:
    tuple val(meta), path('cleaned_counts_filtered.parquet'),                                                         emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    """
    clean_count_data.py \\
        --counts $count_file \\
        --ks-stats $ks_stats_file \\
        --ks-pvalue-threshold $ks_pvalue_threshold
    """

}
