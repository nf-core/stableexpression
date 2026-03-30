process EXPRESSIONATLAS_GETDATA {

    label 'process_single'
    label 'can_fail'

    tag "$accession"

    maxForks 8 // limiting to 8 threads at a time to avoid 429 errors with the Expression Atlas API server

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer']  && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/963bb5cfef2f27d3c5b2a428b18319c65e4d6ff428be08cf3e124e4f9a25a234/data':
        'community.wave.seqera.io/library/bioconductor-expressionatlas_r-base_r-optparse:e15047a6b3701e2c' }"

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
    download_eatlas_data.R --accession $accession
    """

    stub:
    """
    touch acc.raw.counts.csv
    touch acc.design.csv
    """

}
