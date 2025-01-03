process GENE_VARIATION {

    debug true

    publishDir "${params.outdir}/gene_variation"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_files, stageAs: "?/*"
    path metadata_files, stageAs: "?/*"
    path mapping_files, stageAs: "?/*"
    val gene_variation_method

    output:
    path 'gene_variations.csv',                                                                                       emit: csv
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    """
    get_gene_variations.py --counts "$count_files" --metadata "$metadata_files" --mappings "$mapping_files" --method "$gene_variation_method"
    """

}
