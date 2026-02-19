process DOWNLOAD_NCBI_ANNOTATION {

    label 'process_single'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/4709652d855d874806dcd77d35cfe69b0ffec872213cbd573511180f03c096dc/data':
        'community.wave.seqera.io/library/httpx_python_tenacity:2ece9866afa83f4e' }"

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
