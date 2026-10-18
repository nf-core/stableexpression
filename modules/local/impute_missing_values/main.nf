process IMPUTE_MISSING_VALUES {

    label 'process_high'
    label 'process_long'

    tag "${meta.platform}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/79a1048afbb01fd462b7943bdccf16d1ac79fa424a3f9db9f9a76449b5121b88/data':
        'community.wave.seqera.io/library/python_polars_scikit-learn_faiss_pruned:0cb55f1a2dc6ec0e' }"

    input:
    tuple val(meta), path(count_file)

    output:
    tuple val(meta), path('*.imputed.parquet'),                                                                         emit: counts
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                   topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('scikit-learn'), eval('python3 -c "import sklearn; print(sklearn.__version__)"'), topic: versions
    tuple val("${task.process}"), val('faiss'),        eval('python3 -c "import faiss; print(faiss.__version__)"'),     topic: versions
    tuple val("${task.process}"), val('fknni'),        eval('python3 -c "import fknni; print(fknni.__version__)"'),     topic: versions

    script:
    """
    impute_missing_values.py \\
        --counts $count_file \\
        --out ${meta.platform}.imputed.parquet
    """

    stub:
    """
    touch stub.imputed.parquet
    """

}
