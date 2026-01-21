process NORMALISATION_DESEQ2 {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/cef7164b168e74e5db11dcd9acf6172d47ed6753e4814c68f39835d0c6c22f6d/data':
        'community.wave.seqera.io/library/bioconductor-deseq2_r-base_r-optparse:c84cd7ffdb298fa7' }"

    input:
    tuple val(meta), path(count_file), path(design_file)

    output:
    tuple val(meta), path('*.cpm.csv'),                  optional: true,                                             emit: cpm
    tuple val(meta.dataset), path("failure_reason.txt"), optional: true,                                             topic: normalisation_failure_reason
    tuple val(meta.dataset), path("warning_reason.txt"), optional: true,                                             topic: normalisation_warning_reason
    tuple val("${task.process}"), val('R'),      eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
    tuple val("${task.process}"), val('DESeq2'), eval('Rscript -e "cat(as.character(packageVersion(\'DESeq2\')))"'), topic: versions

    script:
    """
    normalise_with_deseq2.R \\
        --counts $count_file \\
        --design $design_file
    """


}
