process EXPRESSIONATLAS_GETACCESSIONS {

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e459ae44332297f0429e7dd501bc3a6f9b5504b13e2db0002a5d3021cc9ac443/data':
        'community.wave.seqera.io/library/nltk_pandas_python_requests_tenacity:a29bfda256e4f39f' }"

    input:
    val species
    val keywords

    output:
    path 'accessions.txt',                                                                                            emit: txt
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions
    tuple val("${task.process}"), val('nltk'),     eval('python3 -c "import nltk; print(nltk.__version__)"'),         topic: versions


    script:
    def keywords_string = keywords.split(',').collect { it.trim() }.join(' ')

    // the folder where nltk will download data needs to be writable (necessary for singularity)
    if (keywords_string == "") {
        """
        NLTK_DATA=$PWD get_eatlas_accessions.py --species $species
        """
    } else {
        """
        NLTK_DATA=$PWD get_eatlas_accessions.py --species $species --keywords $keywords_string
        """
    }


    stub:
    """
    touch accessions.csv
    """

}
