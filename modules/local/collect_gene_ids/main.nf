process COLLECT_GENE_IDS {

    tag "chunk ${task.index}"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/604657081a64b39e17bb6ad307e545aa6aebf4133b64d6766515c9789bb2d304/data':
        'community.wave.seqera.io/library/pandas_tqdm:2ca37c1047243549' }"

    input:
    path count_files, stageAs: "?/*"

    output:
    path 'all_gene_ids.txt',                                                                                          emit: gene_ids
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('tqdm'),     eval('python3 -c "import tqdm; print(tqdm.__version__)"'),         topic: versions

    script:
    """
    collect_gene_ids.py \\
        --counts "$count_files"
    """

}
