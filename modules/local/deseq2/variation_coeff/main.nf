process DESEQ2_VARIATIONCOEFF {

    // label 'error_retry'
    debug true

    conda "${moduleDir}/environment.yml"

    input:
    val count_file


    when:
    task.ext.when == null || task.ext.when

    script:
    template "deseq2.R"


}
