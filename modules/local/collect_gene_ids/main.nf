process COLLECT_GENE_IDS {

    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70c17cde84896904c0620d614cba74ff029f1255db64e66416e63c91b7c959a2/data':
        'community.wave.seqera.io/library/python_tqdm:4e039400f75bdad0' }"

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
