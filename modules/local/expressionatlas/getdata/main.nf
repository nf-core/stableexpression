process EXPRESSIONATLAS_GETDATA {

    // limiting to 8 threads at a time to avoid 429 errors with the Expression Atlas API server
    maxForks 8

    // ignoring accessions that cannot be retrieved from Expression Atlas (the script throws a 100 in this case)
    errorStrategy { task.exitStatus == 100 ? 'ignore' : 'terminate' }

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
