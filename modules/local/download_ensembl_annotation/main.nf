process DOWNLOAD_ENSEMBL_ANNOTATION {

    label 'process_single'

    tag "${species}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/98/980a21a12b628a41a6c08a91d4f6646d1122f0d0e38387f724d4f4ee020b8b1d/data':
        'community.wave.seqera.io/library/bs4_httpx_pandas_python_pruned:13dbe891a99b6884' }"

    input:
    val species

    output:
    path "*.gff3.gz",                                                                                           emit: gff3
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('httpx'),  eval('python3 -c "import httpx; print(httpx.__version__)"'),   topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'), topic: versions
    tuple val("${task.process}"), val('bs4'),    eval('python3 -c "import bs4; print(bs4.__version__)"'),       topic: versions
    tuple val("${task.process}"), val('tqdm'),   eval('python3 -c "import tqdm; print(tqdm.__version__)"'),     topic: versions

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
