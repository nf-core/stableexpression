process DETECT_RARE_GENES {

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path(gene_id_mapping_file)
    path(gene_id_occurrences_file)
    val(nb_datasets)
    val(min_occurrence_frequency)
    val(min_occurrence_quantile)

    output:
    path('valid_gene_ids.txt'),                                                                                       emit: valid_gene_ids
    path('total_gene_id_occurrence_quantiles.csv'),                                                                   topic: total_gene_id_occurrence_quantiles
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    detect_rare_genes.py \\
        --occurrences $gene_id_occurrences_file \\
        --mappings $gene_id_mapping_file \\
        --nb-datasets $nb_datasets \\
        --min-occurrence-frequency $min_occurrence_frequency \\
        --min-occurrence-quantile $min_occurrence_quantile

    """


    stub:
    """
    touch fake.validated_genes.txt
    """

}
