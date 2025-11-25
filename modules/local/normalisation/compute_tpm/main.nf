process NORMALISATION_COMPUTE_TPM {

    label 'process_single'

    tag "${meta.dataset}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3d7126100b0eb7cb53dfb50291707ea8dda3b9738b76551ab73605d0acbe114b/data':
        'community.wave.seqera.io/library/pandas:2.3.3--5a902bf824a79745' }"

    input:
    tuple val(meta), path(count_file)
    path gene_lengths_file

    output:
    tuple val(meta), path('*.tpm.csv'), optional: true,                                                               emit: counts
    tuple val(meta.dataset), path("failure_reason.txt"),    optional: true,                                           topic: normalisation_failure_reason
    tuple val(meta.dataset), path("warning_reason.txt"),    optional: true,                                           topic: normalisation_warning_reason
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions

    script:
    """
    compute_tpm.py \\
        --counts $count_file \\
        --gene-lengths $gene_lengths_file
    """


}
