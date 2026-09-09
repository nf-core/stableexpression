process DETECT_RARE_GENES {

    label 'process_low_requirements'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00f1434368763cebf37466cfaaaf069f971f7eae65b010169975c50d084e5af3/data':
        'community.wave.seqera.io/library/polars_python:1a4a3322c56bfeb9' }"

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
    """
    detect_rare_genes.py \\
        --occurrences $gene_id_occurrences_file \\
        --mappings $gene_id_mapping_file \\
        --nb-datasets $nb_datasets \\
        --min-occurrence-frequency $min_occurrence_frequency \\
        --min-occurrence-quantile $min_occurrence_quantile

    """
}
