process EDGER_NORMALIZE {

    debug true

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.log_cpm.csv'),                                                     emit: csv
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "R.version.string"'),          topic: versions
    tuple val("${task.process}"), val('edgeR'), eval('Rscript -e "packageVersion(\'edgeR\')"'), topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def design_file = meta.design
    """
    edger_normalize.R --counts "$count_file" --design "$design_file"
    """

}
