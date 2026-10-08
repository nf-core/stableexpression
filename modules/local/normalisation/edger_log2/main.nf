process NORMALISATION_EDGER_LOG2 {

    label 'process_single_cpu'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f3/f3f6f1ee98750129d415998fae3312e00c42e09f5e21ab29cca877518aac68cf/data':
        'community.wave.seqera.io/library/bioconductor-edger_r-optparse_r-arrow:2b55e7d643191e09' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.edger_log2.parquet'),            optional: true,                                            emit: counts
    tuple val(meta.dataset), path("failure_reason.txt"), optional: true,                                            topic: normalisation_failure_reason
    tuple val(meta.dataset), path("warning_reason.txt"), optional: true,                                            topic: normalisation_warning_reason
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),  topic: versions
    tuple val("${task.process}"), val('edgeR'), eval('Rscript -e "cat(as.character(packageVersion(\'edgeR\')))"'),  topic: versions

    script:
    """
    edger_log2_normalisation.R \\
        --counts $count_file
    """

    stub:
    """
    touch stub.edger_log2.parquet
    """

}
