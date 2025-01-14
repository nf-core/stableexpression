process EDGER_NORMALISE {

    label 'process_low'

    publishDir "${params.outdir}/normalisation/edger"

    tag "${meta.dataset}"

    // ignoring cases when the count dataframe gets empty after filtering (the script throws a 100 in this case)
    // the subsequent steps will not be run for this dataset
    errorStrategy { task.exitStatus == 100 ? 'ignore' : 'terminate' }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-edger_r-base_r-optparse:aa98b0aff834f649':
        'community.wave.seqera.io/library/bioconductor-edger_r-base_r-optparse:400aaabddeea1574' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.cpm.csv'),                                                                             emit: cpm
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
    tuple val("${task.process}"), val('edgeR'), eval('Rscript -e "cat(as.character(packageVersion(\'edgeR\')))"'),  topic: versions


    script:
    def design_file = meta.design
    """
    edger_normalise.R --counts "$count_file" --design "$design_file"
    """

}
