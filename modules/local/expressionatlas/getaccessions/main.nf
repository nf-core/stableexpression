process EXPRESSIONATLAS_GETACCESSIONS {

    label 'process_low'

    conda "${moduleDir}/spec-file.txt"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f2219a174683388670dc0817da45717014aca444323027480f84aaaf12bfb460/data':
        'community.wave.seqera.io/library/nltk_data_pandas_pyyaml_requests_tenacity:5f5f82f858433879' }"

    input:
    val species
    val keywords

    output:
    path "accessions.txt",                                                                                            emit: accessions
    path "all_experiments.metadata.tsv",                                                                              emit: all_eatlas_experiment_metadata
    path "species_experiments.metadata.tsv",                                                                          topic: species_eatlas_experiment_metadata
    path "filtered_experiments.metadata.tsv", optional: true,                                                         topic: filtered_eatlas_experiment_metadata
    path "filtered_experiments.keywords.yaml", optional: true,                                                        topic: filtered_eatlas_experiment_keywords
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('requests'), eval('python3 -c "import requests; print(requests.__version__)"'), topic: versions
    tuple val("${task.process}"), val('nltk'),     eval('python3 -c "import nltk; print(nltk.__version__)"'),         topic: versions
    tuple val("${task.process}"), val('pyyaml'),   eval('python3 -c "import yaml; print(yaml.__version__)"'),         topic: versions
    tuple val("${task.process}"), val('pandas'),   eval('python3 -c "import pandas; print(pandas.__version__)"'),     topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def keywords_string = keywords.split(',').collect { it.trim() }.join(' ')

    // the folder where nltk will download data needs to be writable (necessary for singularity)
    if (keywords_string == "") {
        """
        NLTK_DATA=$PWD get_eatlas_accessions.py \\
            --species $species
        """
    } else {
        """
        NLTK_DATA=$PWD get_eatlas_accessions.py \\
            --species $species \\
            --keywords $keywords_string
        """
    }

    stub:
    """
    touch accessions.txt all_experiments.metadata.tsv filtered_experiments.metadata.tsv filtered_experiments.keywords.yaml
    """

}
