process GEO_GETDATA {

    label 'process_single'

    tag "$accession"

    maxForks 8 // limiting to 8 threads at a time to avoid 429 errors with the NCBI server

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4cb08d96e62942e7b6288abf2cfd30e813521a022459700e610325a3a7c0b1c8/data':
        'community.wave.seqera.io/library/bioconductor-geoquery_r-base_r-dplyr_r-optparse:fcd002470b7d6809' }"

    input:
    val accession
    val species

    output:
    path("*.counts.csv"),                             optional: true,                                                                   emit: counts
    path("*.design.csv"),                             optional: true,                                                                   emit: design
    path("rejected/**"),                              optional: true,                                                                   emit: rejected
    tuple val(accession), path("failure_reason.txt"), optional: true,                                                                   topic: geo_failure_reason
    tuple val(accession), path("warning_reason.txt"), optional: true,                                                                   topic: geo_warning_reason
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
