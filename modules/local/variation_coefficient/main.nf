process VARIATION_COEFFICIENT {

    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/733cbf61013292f639ec24adcec9548a119ea6254d3ba51c6503ffaba6acda4f/data':
        'community.wave.seqera.io/library/r-base_r-optparse:f7a5d8afb6d6fa3d' }"

    input:
    path count_files

    output:
    path 'variation_coefficients.csv',                                                          emit: csv
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "R.version.string"'),          topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_variation_coefficient.R --count-files "$count_files"
    """

}
