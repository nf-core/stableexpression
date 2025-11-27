process MERGE_COUNTS {

    label "process_high_memory"

    memory { def calc = (dataset_size / 50000).toInteger()
        def result = Math.max(1, calc)  // Ensure at least 1 MB
        def multiplicator = 1 + 0.2 * task.attempt // increase memory usage with each attempt by 20%
        return 1.MB * result * multiplicator
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/90617e987f709570820b8e7752baf9004ba85917111425d4b44b429b27b201ca/data':
        'community.wave.seqera.io/library/polars_tqdm:54b124dde91d1bf3' }"

    input:
    path count_files, stageAs: "?/*"
    val dataset_size

    output:
    path 'all_counts.parquet',                                                                                        emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('tqdm'),     eval('python3 -c "import tqdm; print(tqdm.__version__)"'),         topic: versions

    script:
    """
    merge_counts.py \\
        --counts "$count_files"
    """

}
