process EXPRESSIONATLAS_GETACCESSIONS {

    // label 'error_retry'
    debug true

    afterScript """${baseDir}/bin/write_versions.py ${moduleDir}"""

    conda "${moduleDir}/environment.yml"


    input:
    tuple val(species), val(keywords)

    output:
    path 'accessions.csv'                , emit: accession


    when:
    task.ext.when == null || task.ext.when

    script:
    template "get_expression_atlas_accessions.py"


    stub:
    """
    touch accessions.csv
    touch found.json
    """

}
