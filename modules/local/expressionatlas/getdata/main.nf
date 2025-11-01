process EXPRESSIONATLAS_GETDATA {

    label 'process_single'

    tag "$accession"

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
