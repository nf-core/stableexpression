process VARIATION_COEFFICIENT {

    publishDir "${params.outdir}/variation_coefficients"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ea/eac9296641a6b7b5052cad88b74f09d0f48ae25bca90468645b22bf37bec74b8/data':
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-optparse_r-tibble:afb27df3e35caae2' }"

    input:
    path count_file
    path metadata
    path mapping

    output:
    path 'variation_coefficients.csv',                                                                           emit: csv
    tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),   topic: versions


    script:
    """
    get_variation_coefficients.py --counts "$count_file" --metadata "$metadata" --mapping "$mapping"
    """

}
