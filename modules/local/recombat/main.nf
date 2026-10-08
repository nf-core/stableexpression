process RECOMBAT {

    tag "${meta.platform}"

    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc13b3c613c00dffcab7f37bd14d796e657673b9072def0db1eb3d99f906d6ae/data':
        'community.wave.seqera.io/library/python_polars_scikit-learn:174e3827fae8c638' }"

    input:
    tuple val(meta), path(count_file), path(design_file)

    output:
    tuple val(meta), path('*.corrected.parquet'),                                                                       emit: counts
    tuple val(meta), env('NB_EXCLUDED_ORPHAN_SAMPLES'),                                                                 emit: nb_excluded_orphan_samples
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                   topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval('python3 -c "import sklearn; print(sklearn.__version__)"'), topic: versions

    script:
    def prefix = "${meta.platform}"
    """
    correct_batch_effects_with_recombat.py \\
        --counts $count_file \\
        --design $design_file \\
        --out ${prefix}.corrected.parquet

    NB_EXCLUDED_ORPHAN_SAMPLES=\$(cat nb_excluded_orphan_samples.txt)
    """

    stub:
    """
    touch stub.corrected.parquet
    NB_EXCLUDED_ORPHAN_SAMPLES=0
    """

}
