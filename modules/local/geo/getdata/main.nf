process GEO_GETDATA {

    label 'process_single'

    tag "$accession"

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
