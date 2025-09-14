process GET_CANDIDATE_GENES {

    label 'process_single'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_file
    path stat_file
    val candidate_selection_descriptor
    val nb_top_stable_genes

    output:
    path 'candidate_counts.parquet',                                                                                  emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions

    script:
    """
    get_candidate_genes.py \\
        --counts $count_file \\
        --stats $stat_file \\
        --candidate_selection_descriptor $candidate_selection_descriptor \\
        --nb-top-stable-genes $nb_top_stable_genes
    """

}
