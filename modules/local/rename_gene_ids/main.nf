process RENAME_GENE_IDS {

    label 'process_low'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
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
    tuple val(meta.dataset), env("NB_FINAL"), env("NB_MERGED"), env("NB_UNMAPPED"),                                   topic: id_mapping_stats
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def mapping_arg  = gene_id_mapping_file ? "--mappings $gene_id_mapping_file" : ""
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    rename_gene_ids.py \\
        --count-file "$count_file" \\
        $mapping_arg

    NB_UNMAPPED=\$(cat unmapped.txt)
    NB_MERGED=\$(cat merged.txt)
    NB_FINAL=\$(cat final.txt)
    """


    stub:
    """
    touch fake_renamed.csv
    """

}
