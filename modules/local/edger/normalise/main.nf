process EDGER_NORMALISE {

    publishDir "${params.outdir}/normalisation/edger"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89bbc9544e18b624ed6d0a30e701cf8cec63e063cc9b5243e1efde362fe92228/data':
        'community.wave.seqera.io/library/bioconductor-edger_r-base_r-optparse:400aaabddeea1574' }"

    input:
    tuple val(meta), path(count_file)
    val allow_zero_counts

    output:
    tuple val(meta), path('*.log_cpm.csv'),                                                                         emit: csv
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
    tuple val("${task.process}"), val('edgeR'), eval('Rscript -e "cat(as.character(packageVersion(\'edgeR\')))"'),  topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def design_file = meta.design
    def allow_zeros_arg = allow_zero_counts ? '--allow-zeros' : ''
    """
    edger_normalise.R --counts "$count_file" --design "$design_file" $allow_zeros_arg
    """

}
