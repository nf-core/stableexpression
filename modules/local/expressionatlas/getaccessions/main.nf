process EXPRESSIONATLAS_GETACCESSIONS {

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/nltk_pandas_python_requests_tenacity:3a46f42502407ef2':
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
