process GEO_GETDATA {

    label 'process_single'

    // limiting to 8 threads at a time to avoid 429 errors with the Expression Atlas API server
    maxForks 8

    tag "$accession"

    errorStrategy {
        if (task.exitStatus == 100) {
            // ignoring accessions that cannot be retrieved from GEO
            log.warn("Could not retrieve data for accession ${accession}. This could be a transient network issue or a permission error.")
            return 'ignore'
        } else if (task.exitStatus == 101) {
            log.warn("GEO dataset with accession ${accession} contains multiple files.")
            return 'ignore'
        } else if (task.exitStatus == 110) {
            log.warn("GEO dataset for accession ${accession} does not seem normalised.")
            return 'ignore'
        } else if (task.exitStatus == 111) {
            log.warn("GEO dataset for accession ${accession} seems normalised but not log-transformed.")
            return 'ignore'
        } else if (task.exitStatus == 112) {
            log.warn("GEO dataset for accession ${accession} are of unclear origin. Could not infer normalisation state.")
            return 'ignore'
        } else if (task.exitStatus == 137) { // override default behaviour to sleep some time before retry
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
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4cb08d96e62942e7b6288abf2cfd30e813521a022459700e610325a3a7c0b1c8/data':
        'community.wave.seqera.io/library/bioconductor-geoquery_r-base_r-dplyr_r-optparse:fcd002470b7d6809' }"

    input:
    val accession
    val species

    output:
    path "*.design.csv", optional: true,                                                                                                emit: design
    path "*.counts.csv", optional: true,                                                                                                emit: counts
    tuple val("${task.process}"), val('R'),               eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),            topic: versions
    tuple val("${task.process}"), val('GEOquery'),        eval('Rscript -e "cat(as.character(packageVersion(\'GEOquery\')))"'),         topic: versions
    tuple val("${task.process}"), val('dplyr'),           eval('Rscript -e "cat(as.character(packageVersion(\'dplyr\')))"'),            topic: versions

    script:
    """
    download_geo_data.R \\
        --accession $accession \\
        --species $species
    """

    stub:
    """
    touch acc.microarray.normalised.counts.csv
    touch acc.design.csv
    """

}
