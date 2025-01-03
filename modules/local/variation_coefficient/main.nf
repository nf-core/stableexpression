process VARIATION_COEFFICIENT {

    debug true

    publishDir "${params.outdir}/variation_coefficients"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_files, stageAs: "?/*"
    path metadata_files, stageAs: "?/*"
    path mapping_files, stageAs: "?/*"

    output:
    path 'variation_coefficients.csv',                                                                           emit: csv
    tuple val("${task.process}"), val('R'), eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),   topic: versions


    script:
    """
    get_variation_coefficients.py --counts "$count_files" --metadata "$metadata_files" --mappings "$mapping_files"
    """

}
