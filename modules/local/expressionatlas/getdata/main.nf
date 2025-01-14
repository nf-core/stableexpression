process EXPRESSIONATLAS_GETDATA {

    label 'process_low'

    // limiting to 8 threads at a time to avoid 429 errors with the Expression Atlas API server
    maxForks 8

    tag "$accession"

    errorStrategy = {
        if (task.exitStatus == 100) {
            // ignoring accessions that cannot be retrieved from Expression Atlas (the script throws a 100 in this case)
            return 'ignore'
        } else if (task.exitStatus == 137) { // override default behaviour to sleep some time before retry
            // in case of OOM errors, we wait a bit and try again
            sleep(Math.pow(2, task.attempt) * 2000 as long)
            return 'retry'
        } else {
            return 'terminate'
        }
    }
    maxRetries = 5

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-expressionatlas_r-base_r-optparse:d24070d263d42ce2':
        'community.wave.seqera.io/library/bioconductor-expressionatlas_r-base_r-optparse:ca0f8cd9d3f44af9' }"

    input:
    val(accession)

    output:
    path "*.design.csv",                                                                                         emit: design
    path "*.raw.csv",                   optional: true,                                                          emit: raw
    path "*.normalised.csv",            optional: true,                                                          emit: normalised
    tuple val("${task.process}"), val('R'),               eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),            topic: versions
    tuple val("${task.process}"), val('ExpressionAtlas'), eval('Rscript -e "cat(as.character(packageVersion(\'ExpressionAtlas\')))"'),  topic: versions


    script:
    """
    get_eatlas_data.R --accession $accession
    """

    stub:
    """
    touch acc.raw.csv
    touch acc.design.csv
    """

}
