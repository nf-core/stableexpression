process IDMAPPING {

    // label 'error_retry'
    debug true

    // limiting to 1 thread at a time, otherwise there are two many requests to NCBI
    maxForks 1

    afterScript """${baseDir}/bin/write_versions.py ${moduleDir}"""

    conda "${moduleDir}/environment.yml"

    input:
    tuple path(count_file), val(species)

    output:
    path '*_renamed.csv'                     , emit: csv


    when:
    task.ext.when == null || task.ext.when

    script:
    template "map_ids.py"


    stub:
    """
    touch fake_renamed.csv
    """

}
