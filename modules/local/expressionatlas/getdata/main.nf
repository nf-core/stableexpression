process EXPRESSIONATLAS_GETDATA {

    // when there are network issues, we retry the download with a backoff
    // errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    // maxRetries 5

    // limiting threads to avoid issues with the Expression Atlas API
    maxForks 4

    tag "$accession"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fd21450c3a3f7df37fa0480170780019e9686be319da1c9e10712f7f17cca26/data':
        'community.wave.seqera.io/library/bioconductor-expressionatlas_r-base_r-optparse:ca0f8cd9d3f44af9' }"

    input:
    val(accession)

    output:
    tuple val(accession), path("*.design.csv"), path("*raw.csv"),                   optional: true,                                     emit: raw
    tuple val(accession), path("*.design.csv"), path("*normalised.csv"),            optional: true,                                     emit: normalised
    tuple val("${task.process}"), val('R'),               eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),            topic: versions
    tuple val("${task.process}"), val('ExpressionAtlas'), eval('Rscript -e "cat(as.character(packageVersion(\'ExpressionAtlas\')))"'),  topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_eatlas_data.R --accession $accession
    """

    stub:
    """
    touch acc.raw.csv
    touch acc.design.csv
    """

}
