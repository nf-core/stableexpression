process EXPRESSIONATLAS_GETDATA {

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
    path("*.csv")                             , emit: csv
    path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'get_expression_atlas_data.R'

    stub:
    """
    touch fake.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'R.Version()\$version.string' | sed -n 's|\\[1\\] "R version \\(.*\\) (.*|\\1|p')
        bioconductor-geoquery: \$(Rscript -e 'packageVersion("GEOquery")' | sed -n 's|\\[1\\] ‘\\(.*\\)’|\\1|p')
    END_VERSIONS
    """

}
