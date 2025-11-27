process NORMALISATION_EDGER {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89bbc9544e18b624ed6d0a30e701cf8cec63e063cc9b5243e1efde362fe92228/data':
        'community.wave.seqera.io/library/bioconductor-edger_r-base_r-optparse:400aaabddeea1574' }"

    input:
    tuple val(meta), path(count_file), path(design_file)

    output:
    tuple val(meta), path('*.cpm.csv'),                  optional: true,                                            emit: cpm
    tuple val(meta.dataset), path("failure_reason.txt"), optional: true,                                            topic: normalisation_failure_reason
    tuple val(meta.dataset), path("warning_reason.txt"), optional: true,                                            topic: normalisation_warning_reason
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
    tuple val("${task.process}"), val('edgeR'), eval('Rscript -e "cat(as.character(packageVersion(\'edgeR\')))"'),  topic: versions

    script:
    """
    normalise_with_edger.R \\
        --counts $count_file \\
        --design $design_file
    """

}
