process EDGER_VARIATIONCOEFF {

    debug true

    conda "${moduleDir}/environment.yml"

    input:
    path count_file

    output:
    path '*_normalized.csv', emit: csv


    when:
    task.ext.when == null || task.ext.when

    script:
    template "edger.R"

}
