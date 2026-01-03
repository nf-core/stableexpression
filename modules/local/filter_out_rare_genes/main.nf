process FILTER_OUT_RARE_GENES {

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path(gene_id_mapping_file)
    path(gene_id_occurrences_file)
    val nb_datasets
    val(min_freq_occurrence)

    output:
    path('valid_gene_ids.txt'), optional: true,                                                                     emit: valid_gene_ids
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    get_genes_with_good_occurrence.py \\
        --occurrences $gene_id_occurrences_file \\
        --mappings $gene_id_mapping_file \\
        --nb-datasets $nb_datasets \\
        --min-freq-occurrence $min_freq_occurrence
    """


    stub:
    """
    touch fake.validated_genes.txt
    """

}
