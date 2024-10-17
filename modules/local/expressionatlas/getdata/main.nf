process EXPRESSIONATLAS_GETDATA {

    debug true

    tag "$accession"

    conda "${moduleDir}/environment.yml"

    input:
    val(accession)

    output:
    path '*.design.csv',                                                                                            emit: metadata
    path '*.raw.csv',        optional: true,                                                                        emit: raw
    path '*.normalized.csv', optional: true,                                                                        emit: normalized
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
    touch fake.csv
    """

}
