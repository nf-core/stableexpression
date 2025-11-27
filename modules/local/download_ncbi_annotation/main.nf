process DOWNLOAD_NCBI_ANNOTATION {

    label 'process_single'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c28c8e613c062828aaee4b950029bc90a1a1aa94d5f61016a588c8ec7be8b65/data':
        'community.wave.seqera.io/library/pandas_requests_tenacity:5ba56df089a9d718' }"

    input:
    val species

    output:
    path "*.gff.gz", emit: gff
    tuple val("${task.process}"), val('python'),      eval("python3 --version | sed 's/Python //'"),                          topic: versions
    tuple val("${task.process}"), val('requests'),    eval('python3 -c "import requests; print(requests.__version__)"'),      topic: versions

    script:
    """
    download_latest_ncbi_annotation.py \\
        --species ${species}

    gzip -n *.gff
    """

    stub:
    """
    touch fake.gff3.gz.txt
    """

}
