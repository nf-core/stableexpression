process EXPRESSIONATLAS_GETACCESSIONS {

    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e4/e40fdee15db481a7c9018d85d73fde63235faad794513039a198f3343f2b0e04/data':
        'community.wave.seqera.io/library/nltk_retry_pip_requests:0e24055eb62456ae' }"

    input:
    val species
    val keywords

    output:
    path 'accessions.txt',                                                                                            emit: txt
    tuple val("${task.process}"), val('python'),   eval('python3 --version'),                                         topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions
    tuple val("${task.process}"), val('nltk'),     eval('python3 -c "import nltk; print(nltk.__version__)"'),         topic: versions


    when:
    task.ext.when == null || task.ext.when


    script:

    def keywords_string = keywords.split(',').collect { it.trim() }.join(' ')

    // the folder where nltk will download data needs to be writable (necessary for singularity)
    if (keywords_string == "") {
        """
        NLTK_DATA=$PWD get_expression_atlas_accessions.py --species $species
        """
    } else {
        """
        NLTK_DATA=$PWD get_expression_atlas_accessions.py --species $species --keywords $keywords_string
        """
    }


    stub:
    """
    touch accessions.csv
    """

}
