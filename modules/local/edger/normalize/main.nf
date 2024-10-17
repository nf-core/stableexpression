process EDGER_NORMALIZE {

    debug true

    conda "${moduleDir}/environment.yml"

    input:
    path count_file

    output:
    path '*_normalized.csv',                                                                    emit: csv
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "R.version.string"'),          topic: versions
    tuple val("${task.process}"), val('edgeR'), eval('Rscript -e "packageVersion(\'edgeR\')"'), topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    edger_normalize.R --count-file "$count_file"
    """

}
