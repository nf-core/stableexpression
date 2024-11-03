process DESEQ2_NORMALIZE {

    debug true

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(accession), path(count_file)
    path design_file

    output:
    path '*.log_cpm.csv',                                                                      emit: csv
    tuple val("${task.process}"), val('R'),      eval('Rscript -e "R.version.string"'),           topic: versions
    tuple val("${task.process}"), val('DESeq2'), eval('Rscript -e "packageVersion(\'DESeq2\')"'), topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    deseq2_normalize.R --counts "$count_file" --design "$design_file" --accession "$accession"
    """


}
