process EXPRESSIONATLAS_GETACCESSIONS {

    // label 'error_retry'
    debug true

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-expressionatlas%3A1.30.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-expressionatlas:1.30.0--r43hdfd78af_0' }"

    input:
    tuple val(species), val(keywords)

    output:
    path("accessions.csv")                  , emit: accession_file
    path("found.json")                      , emit: found

    when:
    task.ext.when == null || task.ext.when


    script:

    template 'get_expression_atlas_accessions.py'

    stub:
    """
    touch accessions.csv
    touch found.json

    """

}
