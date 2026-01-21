process DOWNLOAD_GENOME_ANNOTATION {

    label 'process_single'

    tag "$accession"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a6b13690259900baef6865722cb3a319103acc83b5bcab67504c88bde1e3a9f6/data':
        'community.wave.seqera.io/library/ncbi-datasets-cli_unzip:785aabe86637bae4' }"

    input:
    val(accession)

    output:
    path('genomic.gff'), emit: annotation
    tuple val("${task.process}"), val('ncbi-datasets-cli'), eval("datasets --version | sed 's/datasets version: //g'"), topic: versions

    script:
    """
    datasets download genome accession $accession --include gff3

    unzip -o ncbi_dataset.zip
    mv ncbi_dataset/data/${accession}/* .
    """

}
