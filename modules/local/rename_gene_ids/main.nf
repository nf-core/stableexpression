process RENAME_GENE_IDS {

    label 'process_low'

    tag "${meta.dataset}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c28c8e613c062828aaee4b950029bc90a1a1aa94d5f61016a588c8ec7be8b65/data':
        'community.wave.seqera.io/library/pandas_requests_tenacity:5ba56df089a9d718' }"

    input:
    tuple val(meta), path(count_file)
    path gene_id_mapping_file
    path custom_gene_id_mapping_file

    output:
    tuple val(meta), path('*.renamed.csv'),                 optional: true,                                           emit: counts
    tuple val(meta.dataset), path("failure_reason.txt"),    optional: true,                                           topic: renaming_failure_reason
    tuple val(meta.dataset), path("warning_reason.txt"),    optional: true,                                           topic: renaming_warning_reason
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions

    script:
    def mapping_arg  = gene_id_mapping_file ? "--mappings $gene_id_mapping_file" : ""
    def custom_mapping_arg  = custom_gene_id_mapping_file ? "--custom-mappings $custom_gene_id_mapping_file" : ""
    """
    rename_gene_ids.py \\
        --count-file "$count_file" \\
        $mapping_arg \\
        $custom_mapping_arg
    """


    stub:
    """
    touch fake_renamed.csv
    """

}
