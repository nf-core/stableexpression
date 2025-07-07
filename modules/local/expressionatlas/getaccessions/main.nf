process EXPRESSIONATLAS_GETACCESSIONS {

    label 'process_low'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e3d9b407277b8bb8f8850eba40724b1cae9bd6e11ae0019011af82e6ac17cd4/data':
        'community.wave.seqera.io/library/nltk_pandas_python_pyyaml_pruned:2218f9c10723fbf3' }"

    input:
    val species
    val keywords

    output:
    path "accessions.txt",                                                                                            emit: txt
    path "filtered_experiments.yaml",                                                                                 emit: filtered_experiments
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions
    tuple val("${task.process}"), val('nltk'),     eval('python3 -c "import nltk; print(nltk.__version__)"'),         topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def keywords_string = keywords.split(',').collect { it.trim() }.join(' ')

    // the folder where nltk will download data needs to be writable (necessary for singularity)
    if (keywords_string == "") {
        """
        NLTK_DATA=$PWD get_eatlas_accessions.py \
            --species $species \
        """
    } else {
        """
        NLTK_DATA=$PWD get_eatlas_accessions.py \
            --species $species \
            --keywords $keywords_string
        """
    }


    stub:
    """
    touch accessions.txt filtered_experiments.yaml
    """

}
