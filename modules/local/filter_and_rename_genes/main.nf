process FILTER_AND_RENAME_GENES {

    label 'process_low_requirements'

    tag "${meta.dataset}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

    input:
    tuple val(meta), path(count_file)
    path gene_id_mapping_file
    path valid_gene_ids_file

    output:
    tuple val(meta), path('*.renamed.parquet'),             optional: true,                                           emit: counts
    tuple val(meta.dataset), path("failure_reason.txt"),    optional: true,                                           topic: renaming_failure_reason
    tuple val(meta.dataset), path("warning_reason.txt"),    optional: true,                                           topic: renaming_warning_reason
    tuple val(meta.dataset), env("NB_FINAL"), env("NB_MERGED"), env("NB_NOT_VALID"), env("NB_UNMAPPED"),              topic: mqc_id_mapping_stats
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def mapping_arg  = gene_id_mapping_file ? "--mappings $gene_id_mapping_file" : ""
    def valid_ids_arg = valid_gene_ids_file ? "--valid-gene-ids $valid_gene_ids_file" : ""
    """
    filter_and_rename_genes.py \\
        --count-file "$count_file" \\
        $mapping_arg \\
        $valid_ids_arg

    NB_UNMAPPED=\$(cat unmapped.txt)
    NB_MERGED=\$(cat merged.txt)
    NB_NOT_VALID=\$(cat not_valid.txt)
    NB_FINAL=\$(cat final.txt)
    """


    stub:
    """
    touch fake_renamed.csv
    """

}
