process DOWNLOAD_ENSEMBL_ANNOTATION {

    label 'process_single'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5fa11d593e2f2d68c60acc6a00c812793112bff4691754c992fff6b038458604/data':
        'community.wave.seqera.io/library/bs4_pandas_requests_tenacity_tqdm:32f7387852168716' }"

    input:
    val species

    output:
    path "*.gff3.gz", emit: gff3
    tuple val("${task.process}"), val('python'),      eval("python3 --version | sed 's/Python //'"),                          topic: versions
    tuple val("${task.process}"), val('requests'),    eval('python3 -c "import requests; print(requests.__version__)"'),      topic: versions
    tuple val("${task.process}"), val('pandas'),      eval('python3 -c "import pandas; print(pandas.__version__)"'),          topic: versions
    tuple val("${task.process}"), val('bs4'),         eval('python3 -c "import bs4; print(bs4.__version__)"'),                topic: versions
    tuple val("${task.process}"), val('tqdm'),        eval('python3 -c "import tqdm; print(tqdm.__version__)"'),                topic: versions

    script:
    """
    download_latest_ensembl_annotation.py \\
        --species ${species}
    """

    stub:
    """
    touch fake.gff3.gz.txt
    """

}
