process RENAME_GENE_IDS {

    label 'process_low'

    tag "${meta.dataset}"

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9b43e446f2c3b794644fd4c1c86ab09ba0afafc0c02e3fcdf45509ffc89fc4d/data':
        'community.wave.seqera.io/library/pandas_polars:29ea1468b5490a67' }"

    input:
    tuple val(meta), path(count_file)
    path gene_id_mapping_file

    output:
    tuple val(meta), path('*.renamed.csv'),                 optional: true,                                           emit: counts
    tuple val(meta.dataset), path("failure_reason.txt"),    optional: true,                                           topic: renaming_failure_reason
    tuple val(meta.dataset), path("warning_reason.txt"),    optional: true,                                           topic: renaming_warning_reason
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions

    script:
    def mapping_arg  = gene_id_mapping_file ? "--mappings $gene_id_mapping_file" : ""
    """
    rename_gene_ids.py \\
        --count-file "$count_file" \\
        $mapping_arg
    """


    stub:
    """
    touch fake_renamed.csv
    """

}
