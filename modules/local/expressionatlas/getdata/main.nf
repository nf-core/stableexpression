process EXPRESSIONATLAS_GETDATA {

    label 'process_single'

    tag "$accession"

    maxForks 8 // limiting to 8 threads at a time to avoid 429 errors with the Expression Atlas API server

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer']  && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fd21450c3a3f7df37fa0480170780019e9686be319da1c9e10712f7f17cca26/data':
        'community.wave.seqera.io/library/bioconductor-expressionatlas_r-base_r-optparse:ca0f8cd9d3f44af9' }"

    input:
    val accession

    output:
    path("*.counts.csv"),                             optional: true,                                                                   emit: counts
    path("*.design.csv"),                             optional: true,                                                                   emit: design
    tuple val(accession), path("failure_reason.txt"), optional: true,                                                                   topic: eatlas_failure_reason
    tuple val(accession), path("warning_reason.txt"), optional: true,                                                                   topic: eatlas_warning_reason
    tuple val("${task.process}"), val('R'),               eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),            topic: versions
    tuple val("${task.process}"), val('ExpressionAtlas'), eval('Rscript -e "cat(as.character(packageVersion(\'ExpressionAtlas\')))"'),  topic: versions

    script:
    """
    which python
    download_eatlas_data.R --accession $accession
    """

    stub:
    """
    touch acc.raw.counts.csv
    touch acc.design.csv
    """

}
