process COLLECT_GENE_IDS {

    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/80/80143d9f5e0bfe1364e7bf621ca8bb45f707fd48aa1ba3712158fc441d7873b0/data':
        'community.wave.seqera.io/library/tqdm:4.67.1--c1e9fac535191e31' }"

    input:
    path count_files, stageAs: "?/*"

    output:
    path 'unique_gene_ids.txt',                                                                                       emit: unique_gene_ids
    path 'gene_id_occurrences.csv',                                                                                   emit: gene_id_occurrences
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('tqdm'),     eval('python3 -c "import tqdm; print(tqdm.__version__)"'),         topic: versions

    script:
    """
    collect_gene_ids.py \\
        --ids "$count_files"
    """

}
