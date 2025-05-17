process NORMALISATION_EDGER {

    label 'process_low'

    tag "${meta.dataset}"

    // ignoring cases when the count dataframe gets empty after filtering (the script throws a 100 in this case)
    // the subsequent steps will not be run for this dataset
    errorStrategy { task.exitStatus == 100 ? 'ignore' : 'terminate' }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89bbc9544e18b624ed6d0a30e701cf8cec63e063cc9b5243e1efde362fe92228/data':
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
    normalise_with_edger.R --counts "$count_file" --design "$design_file"
    """

}
