process COMPUTE_BASE_STATISTICS {

    label 'process_medium'

    errorStrategy {
        if (task.exitStatus == 100) {
            log.error("No count could be found before merging datasets! Please check the provided accessions and datasets and run again")
            return 'terminate'
        } else if ( task.exitStatus in ((130..145) + 104 + 175) ) { // override default behaviour to sleep some time before retry
            // in case of OOM errors, we wait a bit and try again (2 retries)
            if ( task.attempt <= 2) {
                sleep(Math.pow(2, task.attempt) * 2000 as long)
                return 'retry'
            } else {
                log.error("${accession} caused Out of Memory error multiple times. Ignoring this accession.")
                return 'ignore'
            }
        } else {
            return 'terminate'
        }
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_file
    val platform

    output:
    path '*stats_all_genes.csv',                                                                                      emit: stats
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def args = task.ext.args ?: ''
    if ( platform != [] ) {
        args += " --platform $platform"
    }
    """
    compute_base_statistics.py \\
        --counts $count_file \\
        $args
    """

}
