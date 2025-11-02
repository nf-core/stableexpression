process NORMFINDER   {

    label 'process_high'

    errorStrategy {
        if (task.exitStatus == 100) {
            log.warn("Too few genes to run NormFinder.")
            return 'ignore'
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
            return 'ignore'
        }
    }

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0e/0e0445114887dd260f1632afe116b1e81e02e1acc74a86adca55099469b490d9/data':
        'community.wave.seqera.io/library/numba_numpy_polars_tqdm:6923cfab6fc04dec' }"

    input:
    path count_file
    path design_file

    output:
    path('stability_values.normfinder.csv'),                                                                            emit: stability_values
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                       topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),       topic: versions

    script:
    """
    normfinder.py \
        --counts $count_file \
        --design $design_file
    """

    stub:
    """
    touch stability_values.normfinder.csv
    """

}
