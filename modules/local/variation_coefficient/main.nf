process VARIATION_COEFFICIENT {

    debug true

    conda "${moduleDir}/environment.yml"

    input:
    path count_file

    output:
    path 'variation_coefficients.csv',                                                          emit: csv
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "R.version.string"'),          topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_variation_coefficient.R --counts "$count_file"
    """

}
