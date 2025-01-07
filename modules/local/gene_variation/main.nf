process GENE_VARIATION {

    label 'process_low'

    publishDir "${params.outdir}/gene_variation"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/56c02e34754f3f0a79292e5d4052fe6d2e271714ff9a883097eeef52317cc9d2/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_files, stageAs: "?/*"
    path metadata_files, stageAs: "?/*"
    path mapping_files, stageAs: "?/*"

    output:
    path 'stats_all_genes.csv',                                                                                       emit: stats_all_genes
    path 'stats_most_stable_genes.csv',                                                                               emit: stats_most_stable_genes
    path 'count_summary.csv',                                                                                         emit: count_summary
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    """
    get_gene_variations.py --counts "$count_files" --metadata "$metadata_files" --mappings "$mapping_files"
    """

}
