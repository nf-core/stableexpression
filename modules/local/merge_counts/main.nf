process MERGE_COUNTS {

    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebf9af86c494e5b3e93706ae0d611d28e8788b2e99717f87a1f586ecf30b7067/data':
        'community.wave.seqera.io/library/polars_python_tqdm:ca595df92ae9b061' }"

    input:
    tuple val(meta), path(count_files, stageAs: "?/*")

    output:
    tuple val(meta), path('all_counts.parquet'),                                                                      emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('tqdm'),     eval('python3 -c "import tqdm; print(tqdm.__version__)"'),         topic: versions

    script:
    """
    merge_counts.py \\
        --counts "$count_files"
    """

}
