process EXPRESSIONATLAS_GETDATA {

    // label 'error_retry'
    debug true

    afterScript """${baseDir}/bin/write_versions.py ${moduleDir} --no-python"""

    tag "$accession"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-expressionatlas%3A1.30.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-expressionatlas:1.30.0--r43hdfd78af_0' }"

    input:
    val(accession)

    output:
    path '*.csv'                             , emit: csv


    when:
    task.ext.when == null || task.ext.when

    script:
    template 'get_expression_atlas_data.R'

    stub:
    """
    touch fake.csv
    """

}
