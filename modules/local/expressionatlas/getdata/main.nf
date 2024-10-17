process EXPRESSIONATLAS_GETDATA {

    debug true

    tag "$accession"

    conda "${moduleDir}/environment.yml"

    input:
    val(accession)

    output:
    tuple val(accession), path("*.design.csv"),                                                                     emit: design
    tuple val(accession), path("*.csv"),                                                                            emit: csv
    tuple val("${task.process}"), val('R'),               eval('Rscript -e "R.version.string"'),                    topic: versions
    tuple val("${task.process}"), val('ExpressionAtlas'), eval('Rscript -e "packageVersion(\'ExpressionAtlas\')"'), topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_expression_atlas_data.R --accession $accession
    """

    stub:
    """
    touch acc.raw.csv
    touch acc.normalized.csv
    """

}
