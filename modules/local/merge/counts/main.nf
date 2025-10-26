process MERGE_COUNTS {

    memory { def calc = (dataset_size / 10000).toInteger()
        def result = Math.max(1, calc)  // Ensure at least 1 MB
        def multiplicator = 1 + 0.2 * task.attempt // increase memory usage with each attempt by 20%
        return 1.MB * result * multiplicator
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_files, stageAs: "?/*"
    val dataset_size

    output:
    path 'all_counts.parquet',                                                                                        emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    println task.memory
    """
    merge_counts.py \\
        --counts "$count_files"
    """

}
