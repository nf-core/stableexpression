process EDGER_NORMALIZE {

    debug true

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(accession), path(count_file)
    path design_file

    output:
    tuple val(accession), path('*.log_cpm.csv'),                                                emit: csv
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "R.version.string"'),          topic: versions
    tuple val("${task.process}"), val('edgeR'), eval('Rscript -e "packageVersion(\'edgeR\')"'), topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    edger_normalize.R --counts "$count_file" --design "$design_file" --accession "$accession"
    """

}
