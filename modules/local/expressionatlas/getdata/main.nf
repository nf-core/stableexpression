process EXPRESSIONATLAS_GETDATA {

    label 'process_single'

    // limiting to 8 threads at a time to avoid 429 errors with the Expression Atlas API server
    maxForks 8

    tag "$accession"

    errorStrategy {
        if (task.exitStatus == 100) {
            // ignoring accessions that cannot be retrieved from Expression Atlas (the script throws a 100 in this case)
            // sometimes, some datasets are transiently unavailable from Expression Atlas:
            // we ignore them as there is no point in trying again and again
            // they will be available again soon but we can't know when
            // for some other files, they are simply unavailable for good...
            log.warn("Could not retrieve data for accession ${accession}. This could be a transient network issue or a permission error.")
            return 'ignore'
        } else if (task.exitStatus == 101) {
            // some datasets are not associated with experiment summary
            // we ignore them as there they would be useless for us
            log.warn("Failure to download whole dataset for accession ${accession}. No experiment summary found.")
            return 'ignore'
        } else if (task.exitStatus == 102) {
            // unhandled error: we print an extra message to warn the user
            log.warn("Unhandled error occurred with accession: ${accession}")
            return 'ignore'
        } else if (task.exitStatus == 137) { // override default behaviour to sleep some time before retry
            // in case of OOM errors, we wait a bit and try again (2 retries)
            if ( task.attempt in <= 2) {
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
    container "${ workflow.containerEngine in ['singularity', 'apptainer']  && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fd21450c3a3f7df37fa0480170780019e9686be319da1c9e10712f7f17cca26/data':
        'community.wave.seqera.io/library/bioconductor-expressionatlas_r-base_r-optparse:ca0f8cd9d3f44af9' }"

    input:
    val(accession)

    output:
    path "*.design.csv", optional: true,                                                                                                emit: design
    path "*.counts.csv", optional: true,                                                                                                emit: counts
    tuple val("${task.process}"), val('R'),               eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),            topic: versions
    tuple val("${task.process}"), val('ExpressionAtlas'), eval('Rscript -e "cat(as.character(packageVersion(\'ExpressionAtlas\')))"'),  topic: versions

    script:
    """
    get_eatlas_data.R --accession $accession
    """

    stub:
    """
    touch acc.raw.counts.csv
    touch acc.design.csv
    """

}
